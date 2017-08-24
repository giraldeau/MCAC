

template <int N,class mystorage>
__attribute__((pure)) double& storage_elem<N,mystorage>::operator[](const int i)
{
   return (*Storage)[i][index];
}


/** Default constructor in local storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(void) :
    Storage(new array< vector<double>, N>),
    external_storage(nullptr),
    index(0)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);
}

/** Constructor with external storage */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(mystorage& ext_storage, const int id):
    Storage(ext_storage.Storage),
    external_storage(&ext_storage),
    index(id)
{
}

/** Copy constructor */
template <int N,class mystorage>
storage_elem<N,mystorage>::storage_elem(const storage_elem<N,mystorage>& other):
    Storage(nullptr),
    external_storage(other.external_storage),
    index(other.index)
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
    index(other.index)
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
            (*Storage)[j].erase((*Storage)[j].begin() + index);
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
            (*Storage)[j].erase((*Storage)[j].begin() + index);
        delete Storage;
    }
    Storage = other.Storage;
    external_storage = other.external_storage;
    index = other.index;
    other.Storage = nullptr;
    other.external_storage = nullptr;
    other.index=0;
    return *this;
}






template <int N,class elem>
template<class mylist>
void storage_list<N,elem>::Init(const int _N, mylist& owner)
{
    Destroy();

    Storage = new array< vector<double>, 7>;
    external_storage=nullptr;


    index.assign(_N+1, 0);
    list.assign(_N, nullptr);
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(_N, 0.);

    index[0]=_N;
    const int listSize = _N;
    for (size = 0; size < listSize; size++)
    {
        index[size+1] = size+1;
        list[size] = new elem(owner,size);
    }

}


template <int N,class elem>
 __attribute__((pure)) elem& storage_list<N,elem>::operator[](const int i)
{
    return *list[i-1];
}

template <int N,class elem>
void storage_list<N,elem>::Destroy(void)
{
    if (external_storage==nullptr)
    {
        if (Storage!=nullptr)
        {
            int _N = size;
            for (size = _N; size > 0; size--)
            {
                delete list[size-1];
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
    index(),
    Storage(nullptr),
    external_storage(nullptr),
    size(0)
{}

//storage_list(const int size);

/** Constructor with external storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, int _index[]):
    list(),
    index(),
    Storage(parent.Storage),
    external_storage(&parent),
    size(_index[0])
{
        index.assign(size+1, 0);
        list.assign(size, nullptr);

        index[0]=size;
        const int listSize = size;
        //#pragma omp for simd
        for (int i = 0; i < listSize; i++)
        {
            index[i+1] = _index[i+1];
            int iparent = external_storage->index[index[i+1]]-1;
            list[i] = external_storage->list[iparent];
        }
}

template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, int* _index[], const int start, const int end):
    list(),
    index(),
    Storage(parent.Storage),
    external_storage(&parent),
    size(0)
{

    for(int i=start;i<=end;i++)
        size += _index[i][0];

    index.assign(size+1, 0);
    list.assign(size, nullptr);

    index[0]=size;
    int m=0;
    for(int i=start;i<=end;i++)
    {
    #pragma omp for simd
        for(int j=1;j<=_index[i][0];j++)
        {
            index[m+1] = _index[i][j];
            int iparent = external_storage->index[index[m+1]]-1;
            list[m] = external_storage->list[iparent];
            m++;
        }
    }
}

/** Copy constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other):
    list(other.list),
    index(other.index),
    Storage(nullptr),
    external_storage(other.external_storage),
    size(other.size)
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
    list(other.list),
    index(other.index),
    Storage(other.Storage),
    external_storage(other.external_storage),
    size(other.size)
{
    Storage=nullptr;
    external_storage=nullptr;
    size=0;
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
    index = std::move(other.index);
    size = std::move(other.size);

    other.Storage = nullptr;
    other.external_storage = nullptr;
    other.size=0;

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
    swap(first.size, second.size);
}
