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

public:
    int size;

protected:
    vector < elem* > list;
    vector < int > indexInStorage;

    array< vector<double>, N>* Storage;
    const storage_list* external_storage;

public:

    elem& operator[](const int);

    template<class mylist>
    void Init(const int size, mylist&);

    void merge(storage_list&);
    void remove(elem&);


    /** Default constructor in local storage */
    storage_list(void);
    //storage_list(const int size);

    /** Constructor with external storage */
    storage_list(storage_list& parent,int index[]);
    storage_list(storage_list& parent,int* index[],const int start,const int end);

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

private:
    void Destroy();

};

#include "storagelist.tpp"

#endif // STORAGELIST_H
