

template <int N,class elem>
template<class mylist>
void storage_list<N,elem>::Init(const int _N, mylist& owner)
{
    Destroy();

    Storage = new array< vector<double>, N>;
    external_storage=nullptr;


    indexInStorage.assign(_N+1, 0);
    list.assign(_N, nullptr);
    for (int j=0;j<N;j++)
        (*Storage)[j].assign(_N, 0.);

    indexInStorage[0]=_N;
    const int listSize = _N;
    for (size = 0; size < listSize; size++)
    {
        indexInStorage[size+1] = size+1;
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
    size(0),
    list(),
    indexInStorage(),
    Storage(nullptr),
    external_storage(nullptr)
{}

//storage_list(const int size);

/** Constructor with external storage */
template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, int _index[]):
    size(_index[0]),
    list(),
    indexInStorage(),
    Storage(parent.Storage),
    external_storage(&parent)
{
        indexInStorage.assign(size+1, 0);
        list.assign(size, nullptr);

        indexInStorage[0]=size;
        const int listSize = size;
        //#pragma omp for simd
        for (int i = 0; i < listSize; i++)
        {
            indexInStorage[i+1] = _index[i+1];
            int iparent = external_storage->indexInStorage[indexInStorage[i+1]]-1;
            list[i] = external_storage->list[iparent];
        }
}

template <int N,class elem>
storage_list<N,elem>::storage_list(storage_list<N,elem>& parent, int* _index[], const int start, const int end):
    size(0),
    list(),
    indexInStorage(),
    Storage(parent.Storage),
    external_storage(&parent)
{

    for(int i=start;i<=end;i++)
        size += _index[i][0];

    indexInStorage.assign(size+1, 0);
    list.assign(size, nullptr);

    indexInStorage[0]=size;
    int m=0;
    for(int i=start;i<=end;i++)
    {
    #pragma omp for simd
        for(int j=1;j<=_index[i][0];j++)
        {
            indexInStorage[m+1] = _index[i][j];
            int iparent = external_storage->indexInStorage[indexInStorage[m+1]]-1;
            list[m] = external_storage->list[iparent];
            m++;
        }
    }
}

/** Copy constructor */
template <int N,class elem>
storage_list<N,elem>::storage_list(const storage_list<N,elem>& other):
    size(other.size),
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
    size(other.size),
    list(other.list),
    indexInStorage(other.indexInStorage),
    Storage(other.Storage),
    external_storage(other.external_storage)
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
    indexInStorage = std::move(other.indexInStorage);
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

template <int N,class elem>
void storage_list<N,elem>::merge(storage_list<N,elem>& other)
{
    if(external_storage==nullptr)
    {
        for (int i=0;i<N;i++)
        {
            (*Storage)[i].insert((*Storage)[i].end(),(*other.Storage)[i].begin(),(*other.Storage)[i].end());
        }
        for (int j=1;j<=other.size;j++)
        {
            list.push_back(new elem(*other.list[j]));
            size += 1;
            indexInStorage.push_back(size);
        }
    }
    else
    {
        list.insert(list.end(),other.list.begin(),other.list.end());
        indexInStorage.insert(indexInStorage.end(),other.indexInStorage.begin()+1,other.indexInStorage.end());
        indexInStorage[0]+=other.indexInStorage[0];
        size += other.size;
    }
}

template <int N,class elem>
void storage_list<N,elem>::remove(elem& ToBeRemoved)
{
    const int id = ToBeRemoved.indexInStorage;
    delete list[id];
    list.erase(list.begin()+id);
    indexInStorage.erase(indexInStorage.begin()+id+1);
    for (int i=0;i<N;i++)
    {
        (*Storage)[i].erase((*Storage)[i].begin()+id);
    }
    size--;
}
