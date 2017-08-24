#ifndef STORAGE_H
#define STORAGE_H

#include <iostream>
#include <vector>
#include <array>

using namespace std;

template <int N, class mystorage>
class storage_elem
{
    template<int,class>
    friend class storage_list;

public:

    array< vector<double>, N>* Storage;
    mystorage* external_storage;
    int index;

    double& operator[](const int);

    /** Default constructor in local storage */
    storage_elem(void);

    /** Constructor with external storage */
    storage_elem(mystorage& Storage, const int id);

    /** Copy constructor */
    storage_elem(const storage_elem&);

    /** Move constructor */
    storage_elem (storage_elem&&) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~storage_elem(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    storage_elem& operator= (const storage_elem& other);

    /** Move assignment operator */
    storage_elem& operator= (storage_elem&& other) noexcept;

};


template <int N, class elem>
class storage_list
{
    template<int,class>
    friend class storage_elem;

public:
    vector < elem* > list;
    vector < int > index;

    array< vector<double>, N>* Storage;
    const storage_list* external_storage;

    int size;



    elem& operator[](const int);

    template<class mylist>
    void Init(const int size, mylist&);
    void Destroy();


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
};

#include "storage.tpp"

#endif // STORAGE_H
