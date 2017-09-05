#ifndef STORAGELIST_H
#define STORAGELIST_H

#include <array>
#include <vector>

using namespace std;

template <int N, class elem>
class storage_list
{
    template<int,class>
    friend class storage_elem;

protected:
    int _size;

protected:
    vector < elem* > list;
    vector < int > indexInStorage;

    array< vector<double>, N>* Storage;
    const storage_list* external_storage;

public:

    elem& operator[](const int);
    const elem& operator[](const int) const;


    template<class mylist>
    void Init(const int _size, mylist&);

    void merge(storage_list&);
    void remove(elem&);
    int size() const noexcept;

    /** Default constructor in local storage */
    storage_list(void);
    //storage_list(const int size);

    /** Constructor with external storage */
    storage_list(storage_list& parent,vector<int> index);

    /** Copy constructor */
    storage_list(const storage_list& other);

    /** Move constructor */
    storage_list (storage_list&&) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~storage_list(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    storage_list& operator= (const storage_list& other);

    /** Move assignment operator */
    storage_list& operator= (storage_list&& other) noexcept;

    friend void swap<>(storage_list& first, storage_list& second);

    typename vector<elem*>::iterator begin(void);
    typename vector<elem*>::iterator end(void);
    typename vector<elem*>::const_iterator begin(void) const;
    typename vector<elem*>::const_iterator end(void) const;

private:
    void Destroy();

};



template <int N,class elem>
template<class mylist>
void storage_list<N,elem>::Init(const int _N, mylist& owner)
{
    Destroy();

    Storage = new array< vector<double>, N>;
    external_storage=nullptr;

    //preallocation
    indexInStorage.reserve(_N);
    list.reserve(_N);

    //initialize data
    for (vector<double>& data : (*Storage))
        data.assign(_N, 0.);

    const int listSize = _N;
    for (_size = 0; _size < listSize; _size++)
    {
        indexInStorage.push_back(_size);
        list.push_back(new elem(owner,_size));
    }

}

template <int N,class elem>
int storage_list<N,elem>::size(void) const noexcept
{
    return _size;
}

template <int N,class elem>
 __attribute__((pure)) elem& storage_list<N,elem>::operator[](const int i)
{
    return *list[i];
}

template <int N,class elem>
__attribute__((pure)) const elem& storage_list<N,elem>::operator[](const int i) const
{
 return *list[i];
}

template <int N,class elem>
void storage_list<N,elem>::Destroy(void)
{
    if (external_storage==nullptr)
    {
        if (Storage!=nullptr)
        {
            int _N = _size;
            for (_size = _N-1; _size >= 0; _size--)
            {
                delete list[_size];
            }
            delete Storage;
        }
    }
    else
    {
        Storage = nullptr;
    }
}


/** Default constructor in local storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(void):
    _size(0),
    list(),
    indexInStorage(),
    Storage(nullptr),
    external_storage(nullptr)
{}

/** Constructor with external storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, vector<int> _index):
    _size(_index.size()),
    list(),
    indexInStorage(),
    Storage(parent.Storage),
    external_storage(&parent)
{
        indexInStorage.assign(_size, 0);
        list.assign(_size, nullptr);

        const int listSize = _size;
        //#pragma omp for simd
        for (int i = 0; i < listSize; i++)
        {
            indexInStorage[i] = _index[i];
            int iparent = external_storage->indexInStorage[indexInStorage[i]];
            list[i] = external_storage->list[iparent];
        }
}

/** Copy constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other):
    _size(other._size),
    list(other.list),
    indexInStorage(other.indexInStorage),
    Storage(nullptr),
    external_storage(other.external_storage)
{
    if(external_storage==nullptr)
    {
        Storage = new array< vector<double>, N>;
        for (int i=0;i<N;i++)
            Storage[i]=other.Storage[i];
    }
    else
        Storage = other.Storage;
}

/** Move constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list (storage_list&& other) noexcept :/* noexcept needed to enable optimizations in containers */
    _size(other._size),
    list(other.list),
    indexInStorage(other.indexInStorage),
    Storage(other.Storage),
    external_storage(other.external_storage)
{
    Storage=nullptr;
    external_storage=nullptr;
    _size=0;
}

/** Destructor */
template <int N,class elem>
storage_list<N,elem>::~storage_list(void) noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    Destroy();
}

/** Copy assignment operator */
template <int N,class elem>
storage_list<N,elem>& storage_list<N,elem>::operator= (const storage_list<N,elem>& other)
{
    storage_list<N,elem> tmp(other);         // re-use copy-constructor
    *this = std::move(tmp);                  // re-use move-assignment
    return *this;
}

/** Move assignment operator */
template <int N,class elem>
storage_list<N,elem>& storage_list<N,elem>::operator= (storage_list<N,elem>&& other) noexcept
{
    Destroy();

    Storage = other.Storage;
    external_storage = other.external_storage;

    list = std::move(other.list);
    indexInStorage = std::move(other.indexInStorage);
    _size = std::move(other._size);

    other.Storage = nullptr;
    other.external_storage = nullptr;
    other._size=0;

    return *this;
}


template <int N,class elem>
void swap(storage_list<N,elem>& first, storage_list<N,elem>& second)
{
    using std::swap;
    swap(first.list, second.list);
    swap(first.index, second.index);
    swap(first.Storage, second.Storage);
    swap(first.external_storage, second.external_storage);
    swap(first._size, second._size);
}

template <int N,class elem>
void storage_list<N,elem>::merge(storage_list<N,elem>& other)
{
    if(external_storage==nullptr)
    {
        for (int i=0;i<N;i++)
        {
            (*Storage)[i].insert((*Storage)[i].end(),(*other.Storage)[i].begin(),(*other.Storage)[i].end());
        }
        for (const auto& data : other.list)
        {
            list.push_back(new elem(*data));
            indexInStorage.push_back(_size);
            _size += 1;
        }
    }
    else
    {
        list.insert(list.end(),other.list.begin(),other.list.end());
        indexInStorage.insert(indexInStorage.end(),other.indexInStorage.begin(),other.indexInStorage.end());
        _size += other._size;
    }
}

template <int N,class elem>
void storage_list<N,elem>::remove(elem& ToBeRemoved)
{
    const int id = ToBeRemoved.indexInStorage;
    delete list[id];
    list.erase(list.begin()+id);
    indexInStorage.erase(indexInStorage.begin()+id);
    for (int i=0;i<N;i++)
    {
        (*Storage)[i].erase((*Storage)[i].begin()+id);
    }
    _size--;
}


template <int N,class elem>
typename vector<elem*>::iterator storage_list<N,elem>::begin(void)
{
    return list.begin();
}

template <int N,class elem>
typename vector<elem*>::iterator storage_list<N,elem>::end(void)
{
    return list.end();
}

template <int N,class elem>
typename vector<elem*>::const_iterator storage_list<N,elem>::begin(void) const
{
    return list.begin();
}

template <int N,class elem>
typename vector<elem*>::const_iterator storage_list<N,elem>::end(void) const
{
    return list.end();
}

#endif // STORAGELIST_H
