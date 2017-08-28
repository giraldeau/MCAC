

template <int N,class mystorage>
__attribute__((pure)) double& storage_elem<N,mystorage>::operator[](const int i)
{
   return (*Storage)[i][indexInStorage];
}


/** Default constructor in local storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(void) :
    Storage(new array< vector<double>, N>),
    external_storage(nullptr),
    indexInStorage(0)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);
}

/** Constructor with external storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(mystorage& ext_storage, const int id):
    Storage(ext_storage.Storage),
    external_storage(&ext_storage),
    indexInStorage(id)
{
}

/** Copy constructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(const storage_elem<N,mystorage>& other):
    Storage(nullptr),
    external_storage(other.external_storage),
    indexInStorage(other.indexInStorage)
{
    if(external_storage==nullptr)
        Storage = new array< vector<double>, N>;
    else
        Storage = other.Storage;
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
    other.index = 0;

}

/** Destructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::~storage_elem(void) noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(external_storage==nullptr && Storage!=nullptr)
    {
        for (int j=0;j<=6;j++)
            (*Storage)[j].erase((*Storage)[j].begin() + indexInStorage);
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
        for (int j=0;j<=6;j++)
            (*Storage)[j].erase((*Storage)[j].begin() + indexInStorage);
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
