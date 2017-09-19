#ifndef STORAGELIST_H
#define STORAGELIST_H

#include <array>
#include <vector>

template <int N, class elem>
class storage_list
{
    template<int,class>
    friend class storage_elem;

protected:
    std::vector < elem* > list;

    std::array< std::vector<double>, N>* Storage;
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

    /** Constructor with external storage */
    storage_list(storage_list& parent,std::vector<int> index);

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

    friend void std::swap<>(storage_list& first, storage_list& second);

    typename std::vector<elem*>::iterator begin(void);
    typename std::vector<elem*>::iterator end(void);
    typename std::vector<elem*>::const_iterator begin(void) const;
    typename std::vector<elem*>::const_iterator end(void) const;

private:
    void Destroy();

};



template <int N,class elem>
template<class mylist>
void storage_list<N,elem>::Init(const int _N, mylist& owner)
{
    Destroy();

    Storage = new std::array< std::vector<double>, N>;
    external_storage=nullptr;

    //preallocation
    list.reserve(_N);

    //initialize data
    for (std::vector<double>& data : (*Storage))
        data.assign(_N, 0.);

    const int listSize = _N;
    for (int _size = 0; _size < listSize; _size++)
    {
        list.push_back(new elem(owner,_size));
    }

}

template <int N,class elem>
int storage_list<N,elem>::size(void) const noexcept
{
    return int(list.size());
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
            for (int _N = size()-1; _N >= 0; _N--)
            {
                delete list[_N];
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
    list(),
    Storage(nullptr),
    external_storage(nullptr)
{}

/** Constructor with external storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, std::vector<int> _index):
    list(),
    Storage(parent.Storage),
    external_storage(&parent)
{
        list.assign(_index.size(), nullptr);

        const int listSize = size();
        //#pragma omp for simd
        for (int i = 0; i < listSize; i++)
        {
            list[i] = external_storage->list[_index[i]];
        }
}

/** Copy constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other):
    list(),
    Storage(new std::array< std::vector<double>, N>),
    external_storage(nullptr)
{
    for (int i=0;i<N;i++)
        (*Storage)[i].assign((*other.Storage)[i].begin(),(*other.Storage)[i].end());
    list.reserve(other.size());
    const int listSize = other.size();
    for (int _size = 0; _size < listSize; _size++)
    {
        list.push_back(new elem(*other.list[_size]));
    }
}

/** Move constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list (storage_list&& other) noexcept :/* noexcept needed to enable optimizations in containers */
    list(other.list),
    Storage(other.Storage),
    external_storage(other.external_storage)
{
    Storage=nullptr;
    external_storage=nullptr;
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
    *this = move(tmp);                  // re-use move-assignment
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

    other.Storage = nullptr;
    other.external_storage = nullptr;

    return *this;
}


template <int N,class elem>
void swap(storage_list<N,elem>& first, storage_list<N,elem>& second)
{
    using std::swap;
    swap(first.list, second.list);
    swap(first.Storage, second.Storage);
    swap(first.external_storage, second.external_storage);
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
        }
    }
    else
    {
        list.insert(list.end(),other.list.begin(),other.list.end());
    }
}

template <int N,class elem>
void storage_list<N,elem>::remove(elem& ToBeRemoved)
{
    const int id = ToBeRemoved.indexInStorage;
    delete list[id];
    list.erase(list.begin()+id);
    for (int i=0;i<N;i++)
    {
        (*Storage)[i].erase((*Storage)[i].begin()+id);
    }
}


template <int N,class elem>
typename std::vector<elem*>::iterator storage_list<N,elem>::begin(void)
{
    return list.begin();
}

template <int N,class elem>
typename std::vector<elem*>::iterator storage_list<N,elem>::end(void)
{
    return list.end();
}

template <int N,class elem>
typename std::vector<elem*>::const_iterator storage_list<N,elem>::begin(void) const
{
    return list.begin();
}

template <int N,class elem>
typename std::vector<elem*>::const_iterator storage_list<N,elem>::end(void) const
{
    return list.end();
}

#endif // STORAGELIST_H
