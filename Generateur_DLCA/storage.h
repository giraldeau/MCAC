#ifndef STORAGE_H
#define STORAGE_H

#include <array>
#include <vector>

using namespace std;

template <int N, class mystorage>
class storage_elem
{
    template<int,class>
    friend class storage_list;

protected:

    array< vector<double>, N>* Storage;
    mystorage* external_storage;
    int indexInStorage;

public:

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

#include "storage.tpp"

#endif // STORAGE_H
