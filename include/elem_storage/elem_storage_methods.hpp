#ifndef INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP
#define INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP 1
#include "elem_storage/elem_storage.hpp"
#include <array>
#include <vector>


namespace MCAC {
template<int N, class mystorage>
void ElemStorage<N, mystorage>::decrease_index() noexcept {
    index_in_storage--;
}
/** Default constructor in local storage */
template<int N, class mystorage>
ElemStorage<N, mystorage>::ElemStorage() noexcept:
    storage(new std::array<std::vector<double>, N>),
    external_storage(nullptr),
    index_in_storage(0) {
    for (size_t j = 0; j < N; j++) {
        (*storage)[j].assign(1, 0.);
    }
}
/** Constructor with external storage */
template<int N, class mystorage>
ElemStorage<N, mystorage>::ElemStorage(mystorage &ext_storage, const size_t id) noexcept:
    storage(ext_storage.storage),
    external_storage(&ext_storage),
    index_in_storage(id) {
}
/** Destructor */
template<int N, class mystorage>
ElemStorage<N, mystorage>::~ElemStorage() noexcept {
    if (!external_storage && storage) {
        for (size_t j = 0; j < N; j++) {
            (*storage)[j].erase((*storage)[j].begin() + long(index_in_storage));
        }
        delete storage;
    }
}
/** Copy constructor */
template<int N, class mystorage>
template<class elem>
ElemStorage<N, mystorage>::ElemStorage(const ElemStorage<N, mystorage> &other,
                                       elem &sphere,
                                       mystorage &_Storage) noexcept:
    storage(_Storage.storage),
    external_storage(&_Storage),
    index_in_storage(external_storage->size()) {
    for (size_t j = 0; j < N; j++) {
        (*storage)[j].push_back((*other.storage)[j][other.index_in_storage]);
    }
    external_storage->list.push_back(&sphere);
}
///** Copy constructor */
//template<int N, class mystorage>
//ElemStorage<N, mystorage>::ElemStorage(const ElemStorage<N, mystorage> &other) noexcept:
//    storage(new std::array<std::vector<double>, N>),
//    external_storage(nullptr),
//    index_in_storage(0) {
//    for (size_t j = 0; j < N; j++) {
//        (*storage)[j].assign(1, (*other.storage)[j][other.index_in_storage]);
//    }
//}
///** Move constructor */
//template<int N, class mystorage>
//ElemStorage<N, mystorage>::ElemStorage(ElemStorage<N, mystorage> &&other) noexcept :
//    storage(other.storage),
//    external_storage(other.external_storage),
//    index_in_storage(other.index_in_storage) {
//    other.storage = nullptr;
//    other.external_storage = nullptr;
//    other.index_in_storage = 0;
//}
///** Copy assignment operator */
//template<int N, class mystorage>
//ElemStorage<N, mystorage> &ElemStorage<N, mystorage>::operator=(const ElemStorage<N, mystorage> &other) noexcept {
//    ElemStorage<N, mystorage> tmp(other);         // re-use copy-constructor
//    *this = std::move(tmp);                       // re-use move-assignment
//    return *this;
//}
///** Move assignment operator */
//template<int N, class mystorage>
//ElemStorage<N, mystorage> &ElemStorage<N, mystorage>::operator=(ElemStorage<N, mystorage> &&other) noexcept {
//    if (!external_storage && storage) {
//        for (size_t j = 0; j < N; j++) {
//            (*storage)[j].erase((*storage)[j].begin() + long(index_in_storage));
//        }
//        delete storage;
//    }
//    storage = other.storage;
//    external_storage = other.external_storage;
//    index_in_storage = other.index_in_storage;
//    other.storage = nullptr;
//    other.external_storage = nullptr;
//    other.index_in_storage = 0;
//    return *this;
//}
}  // namespace MCAC
#endif // INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP_
