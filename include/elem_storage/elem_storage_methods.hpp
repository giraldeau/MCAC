/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP
#define INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP 1
#include "elem_storage/elem_storage.hpp"
#include <array>
#include <vector>


namespace mcac {
template<int N, class mystorage>
void ElemStorage<N, mystorage>::decrease_index() noexcept {
    index_in_storage--;
}
/** Default constructor in local storage */
template<int N, class mystorage>
ElemStorage<N, mystorage>::ElemStorage() noexcept:
    storage(std::make_shared<std::array<std::vector<double>, N>>()),
    external_storage(nullptr),
    index_in_storage(0) {
    for (size_t j = 0; j < N; j++) {
        (*storage)[j].assign(1, 0.);
    }
}
/** Constructor with external storage */
template<int N, class mystorage>
ElemStorage<N, mystorage>::ElemStorage(mystorage &ext_storage, size_t id) noexcept:
    storage(ext_storage.storage),
    external_storage(&ext_storage),
    index_in_storage(id) {
}
/** Destructor */
template<int N, class mystorage>
ElemStorage<N, mystorage>::~ElemStorage() noexcept {
//    if (!external_storage && storage) {
//        for (size_t j = 0; j < N; j++) {
//            (*storage)[j].erase((*storage)[j].begin() + long(index_in_storage));
//        }
//        delete storage;
//    }
}
/** Copy constructor */
template<int N, class mystorage>
template<class elem>
ElemStorage<N, mystorage>::ElemStorage(const ElemStorage<N, mystorage> &other,
                                       elem &sphere,
                                       mystorage &ext_storage) noexcept:
    storage(ext_storage.storage),
    external_storage(&ext_storage),
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
}  // namespace mcac
#endif //INCLUDE_ELEM_STORAGE_ELEM_STORAGE_METHODS_HPP
