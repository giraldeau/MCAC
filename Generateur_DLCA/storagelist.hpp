#ifndef STORAGELIST_H
#define STORAGELIST_H 1

#define UNUSED(expr) do { (void)(expr); } while (0)

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

    elem& operator[](size_t i);
    const elem& operator[](size_t i) const;


    template<class mylist>
    void Init(size_t _size, mylist& owner);

    void merge(storage_list& other);
    void remove(const size_t& id);
    size_t size() const noexcept;

    template <class mylist>
    elem* add(const elem& other, mylist& owner);


    /** Default constructor in local storage */
    storage_list();

    /** Constructor with external storage */
    explicit storage_list(storage_list& parent,std::vector<size_t> _index);

    /** Copy constructor with external storage */
    template<class mylist>
    storage_list(const storage_list& other, mylist& owner);

    template<class mylist>
    storage_list(const storage_list& other, mylist& owner, mylist& _Storage);


    /** Copy constructor in local storage */
    storage_list(const storage_list& other);

    /** Move constructor */
    storage_list (storage_list&& other) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~storage_list() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    storage_list& operator= (const storage_list& other);

    /** Move assignment operator */
    storage_list& operator= (storage_list&& other) noexcept;

//    friend void std::swap<>(storage_list& first, storage_list& second);

    typename std::vector<elem*>::iterator begin();
    typename std::vector<elem*>::iterator end();
    typename std::vector<elem*>::const_iterator begin() const;
    typename std::vector<elem*>::const_iterator end() const;

private:
    void Destroy();

};



template <int N,class elem>
template<class mylist>
void storage_list<N,elem>::Init(const size_t _size, mylist& owner)
{
    Destroy();

    Storage = new std::array< std::vector<double>, N>;
    external_storage=nullptr;

    //preallocation
    list.reserve(_size);

    //initialize data
    for (std::vector<double>& data : (*Storage))
    {
        data.assign(_size, 0.);
    }

    const size_t listSize = _size;
    for (size_t i = 0; i < listSize; i++)
    {
        list.push_back(new elem(owner,i));
    }

}

template <int N,class elem>
__attribute__((pure)) size_t storage_list<N,elem>::size() const noexcept
{
    return list.size();
}

template <int N,class elem>
 __attribute__((pure)) elem& storage_list<N,elem>::operator[](const size_t i)
{
    return *list[i];
}

template <int N,class elem>
__attribute__((pure)) const elem& storage_list<N,elem>::operator[](const size_t i) const
{
 return *list[i];
}

template <int N,class elem>
void storage_list<N,elem>::Destroy()
{
    if (external_storage==nullptr)
    {
        if (Storage!=nullptr)
        {
            if (size()>0)
            {
                for (size_t _N = size(); _N --> 0;)
                {
                    delete list[_N];
                }
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
storage_list<N,elem>::storage_list():
    list(),
    Storage(nullptr),
    external_storage(nullptr)
{}

/** Constructor with external storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, std::vector<size_t> _index):
    list(),
    Storage(parent.Storage),
    external_storage(&parent)
{
        list.assign(_index.size(), nullptr);

        const size_t listSize = size();
        //#pragma omp for simd
        for (size_t i = 0; i < listSize; i++)
        {
            list[i] = external_storage->list[_index[i]];
        }
}

/** Copy constructor */
template <int N,class elem>
template<class mylist>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other, mylist& owner):
    list(),
    Storage(new std::array< std::vector<double>, N>),
    external_storage(nullptr)
{
    for (size_t i=0;i<N;i++)
    {
        (*Storage)[i].reserve(other.size());
    }
    const size_t listSize = other.size();
    for (size_t _size = 0; _size < listSize; _size++)
    {
        for (size_t i=0;i<N;i++)
        {
            (*Storage)[i][_size] = (*other.Storage)[i][other.list[_size]->indexInStorage];
        }
        list.push_back(new elem(owner,_size));
    }
}
/** Copy constructor */
template <int N,class elem>
template<class mylist>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other, mylist& owner, mylist& _Storage):
    list(),
    Storage(_Storage.Storage),
    external_storage(&_Storage)
{
    UNUSED(owner);

    const size_t start = (*Storage)[0].size();
    for (size_t i=0;i<N;i++)
    {
        (*Storage)[i].reserve(other.size()+start);
    }
    list.reserve(other.size()+start);
    const size_t listSize = other.size();
    for (size_t _size = 0; _size < listSize; _size++)
    {
        list.push_back(new elem(*other.list[_size],_Storage,_size+start));
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
storage_list<N,elem>::~storage_list() noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
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
        for (size_t i=0;i<N;i++)
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
void storage_list<N,elem>::remove(const size_t& id)
{
    delete list[id];
    list.erase(list.begin()+long(id));
    for (size_t i=0;i<N;i++)
    {
        (*Storage)[i].erase((*Storage)[i].begin()+long(id));
    }
    for (size_t i = id; i<list.size();i++)
    {
        list[i]->DecreaseLabel();
    }
}


template <int N,class elem>
typename std::vector<elem*>::iterator storage_list<N,elem>::begin()
{
    return list.begin();
}

template <int N,class elem>
typename std::vector<elem*>::iterator storage_list<N,elem>::end()
{
    return list.end();
}

template <int N,class elem>
typename std::vector<elem*>::const_iterator storage_list<N,elem>::begin() const
{
    return list.begin();
}

template <int N,class elem>
typename std::vector<elem*>::const_iterator storage_list<N,elem>::end() const
{
    return list.end();
}


template <int N,class elem>
template <class mylist>
elem* storage_list<N,elem>::add(const elem& other, mylist& owner)
{
    return new elem(other,owner);
}
#endif // STORAGELIST_H
