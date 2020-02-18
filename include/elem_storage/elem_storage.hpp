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
#ifndef INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP
#define INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP 1
#include <array>
#include <vector>


namespace mcac {
template<int N, class mystorage>
class ElemStorage {
    template<int, class>
    friend
    class StorageList;

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
    ElemStorage(const ElemStorage &other, elem &sphere, mystorage &storage) noexcept;
    explicit ElemStorage(const ElemStorage &other) noexcept = delete;
    /** Move constructor */
    explicit ElemStorage(ElemStorage &&other) noexcept = delete;
    /** Copy assignment operator */
    ElemStorage &operator=(const ElemStorage &other) noexcept = delete;
    /** Move assignment operator */
    ElemStorage &operator=(ElemStorage &&other) noexcept = delete;
};
}  // namespace mcac

#include "elem_storage/elem_storage_methods.hpp"


#endif //INCLUDE_ELEM_STORAGE_ELEM_STORAGE_HPP
