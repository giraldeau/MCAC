#ifndef STORAGE_H
#define STORAGE_H

#include <array>
#include <vector>

template <int N, class mystorage>
class storage_elem
{
    template<int,class>
    friend class storage_list;

protected:

    std::array< std::vector<double>, N>* Storage;
    mystorage* external_storage;
    size_t indexInStorage{0};

public:

    void DecreaseIndex();

    /** Default constructor in local storage */
    storage_elem();

    /** Constructor with external storage */
    explicit storage_elem(mystorage& ext_storage, size_t id);

    /** Copy constructor */
    storage_elem(const storage_elem& other);

    template <class elem>
    storage_elem(const storage_elem& other,elem& sphere, mystorage& _Storage);

    /** Move constructor */
    storage_elem (storage_elem&& other) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~storage_elem() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    storage_elem& operator= (const storage_elem& other);

    /** Move assignment operator */
    storage_elem& operator= (storage_elem&& other) noexcept;

};


template <int N,class mystorage>
void storage_elem<N,mystorage>::DecreaseIndex()
{
    indexInStorage--;
}


/** Default constructor in local storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem() :
    Storage(new std::array< std::vector<double>, N>),
    external_storage(nullptr),
    indexInStorage(0)
{
    for (size_t j=0;j<N;j++)
    {
        (*Storage)[j].assign(1, 0.);
    }
}

/** Constructor with external storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(mystorage& ext_storage, const size_t id):
    Storage(ext_storage.Storage),
    external_storage(&ext_storage),
    indexInStorage(id)
{
}

/** Copy constructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(const storage_elem<N,mystorage>& other):
    Storage(new std::array< std::vector<double>, N>),
    external_storage(nullptr),
    indexInStorage(0)
{
    for (size_t j=0;j<N;j++)
    {
        (*Storage)[j].assign(1, (*other.Storage)[j][other.indexInStorage]);
    }
}
/** Copy constructor */
template <int N,class mystorage>
template <class elem>
storage_elem<N,mystorage>::storage_elem(const storage_elem<N,mystorage>& other, elem& sphere, mystorage& _Storage):
    Storage(_Storage.Storage),
    external_storage(&_Storage),
    indexInStorage(external_storage->size())
{
    for (size_t j=0;j<N;j++)
    {
        (*Storage)[j].push_back((*other.Storage)[j][other.indexInStorage]);
    }
    external_storage->list.push_back(&sphere);
}


/** Move constructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem (storage_elem<N,mystorage>&& other) noexcept :/* noexcept needed to enable optimizations in containers */
    Storage(other.Storage),
    external_storage(other.external_storage),
    indexInStorage(other.indexInStorage)
{
    other.Storage = nullptr;
    other.external_storage = nullptr;
    other.indexInStorage = 0;
}

/** Destructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::~storage_elem() noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(external_storage==nullptr && Storage!=nullptr)
    {
        for (size_t j=0;j<N;j++)
        {
            (*Storage)[j].erase((*Storage)[j].begin() + long(indexInStorage));
        }
        delete Storage;
    }
}



/** Copy assignment operator */
template <int N,class mystorage>
storage_elem<N,mystorage>& storage_elem<N,mystorage>::operator=(const storage_elem<N,mystorage>& other)
{
    storage_elem<N,mystorage> tmp(other);         // re-use copy-constructor
    *this = std::move(tmp);                       // re-use move-assignment
    return *this;
}

/** Move assignment operator */
template <int N,class mystorage>
storage_elem<N,mystorage>& storage_elem<N,mystorage>::operator=(storage_elem<N,mystorage>&& other) noexcept
{
    if(external_storage==nullptr && Storage!=nullptr)
    {
        for (size_t j=0;j<N;j++)
        {
            (*Storage)[j].erase((*Storage)[j].begin() + long(indexInStorage));
        }
        delete Storage;
    }
    Storage = other.Storage;
    external_storage = other.external_storage;
    indexInStorage = other.indexInStorage;
    other.Storage = nullptr;
    other.external_storage = nullptr;
    other.indexInStorage=0;
    return *this;
}

#endif // STORAGE_H
