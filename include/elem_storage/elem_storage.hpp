#ifndef INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP_
#define INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP_ 1
#include <array>
#include <vector>


namespace MCAC {
template<int N, class mystorage>
class ElemStorage {
    template<int, class>
    friend
    class storage_list;

protected:
    std::array<std::vector<double>, N> *storage;
    mystorage *external_storage;
    size_t index_in_storage{0};
public:
    void decrease_index() noexcept;
    /** Default constructor in local storage */
    ElemStorage() noexcept;
    /** Constructor with external storage */
    ElemStorage(mystorage &ext_storage, size_t id) noexcept;
    /** Destructor */
    ~ElemStorage() noexcept;
    /** Copy constructor */
    template<class elem>
    ElemStorage(const ElemStorage &other, elem &sphere, mystorage &_Storage) noexcept;
    explicit ElemStorage(const ElemStorage &other) noexcept = delete;
    /** Move constructor */
    explicit ElemStorage(ElemStorage &&other) noexcept = delete;
    /** Copy assignment operator */
    ElemStorage &operator=(const ElemStorage &other) noexcept = delete;
    /** Move assignment operator */
    ElemStorage &operator=(ElemStorage &&other) noexcept = delete;
};
}  // namespace MCAC

#include "elem_storage/elem_storage_methods.hpp"


#endif //INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP_
