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
#ifndef INCLUDE_LIST_STORAGE_LIST_STORAGE_HPP
#define INCLUDE_LIST_STORAGE_LIST_STORAGE_HPP 1
#include <array>
#include <vector>


namespace mcac {
template<int N, class elem>
class ListStorage {
    template<int, class>
    friend
    class ElemStorage;

protected:
    std::vector<std::shared_ptr<elem>> list;
    std::shared_ptr<std::array<std::vector<double>, N>> storage;
    const ListStorage *external_storage;
public:
//    elem &operator[](size_t i) noexcept;
//    const elem &operator[](size_t i) const noexcept;
    std::shared_ptr<elem> &operator[](size_t i) noexcept;
    const std::shared_ptr<elem> &operator[](size_t i) const noexcept;
    template<class mylist>
    void init(size_t size, mylist &owner) noexcept;
    void merge(ListStorage &other) noexcept;
    void remove(const size_t &id) noexcept;
    size_t size() const noexcept;
    template<class mylist>
    void add(size_t, mylist &owner) noexcept;
    template<class mylist>
    std::shared_ptr<elem> add(const elem &other, mylist &owner) noexcept;
    typename std::vector<std::shared_ptr<elem>>::iterator begin() noexcept;
    typename std::vector<std::shared_ptr<elem>>::iterator end() noexcept;
    typename std::vector<std::shared_ptr<elem>>::const_iterator begin() const noexcept;
    typename std::vector<std::shared_ptr<elem>>::const_iterator end() const noexcept;
private:
    void destroy() noexcept;
public:
    /** Default constructor in local storage */
    ListStorage() noexcept;
    /** Constructor with external storage */
    ListStorage(ListStorage &parent, const std::vector<size_t> &index) noexcept;
    /** Destructor */
    ~ListStorage() noexcept;
    /** Copy constructor with external storage */
    template<class mylist>
    ListStorage(const ListStorage &other, mylist &storage) noexcept;
    /** Move constructor */
    ListStorage(ListStorage &&other) noexcept;
    /** Move assignment operator */
    ListStorage &operator=(ListStorage &&other) noexcept;
    /** Copy constructor in local storage */
    ListStorage(const ListStorage &other) noexcept = delete;
    /** Copy assignment operator */
    ListStorage &operator=(const ListStorage &other) noexcept = delete;
};
}  // namespace mcac
#include "list_storage/list_storage_methods.hpp"


#endif //INCLUDE_LIST_STORAGE_LIST_STORAGE_HPP
