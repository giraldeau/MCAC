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
    std::vector<elem *> list;
    std::array<std::vector<double>, N> *storage;
    const ListStorage *external_storage;
public:
    elem &operator[](size_t i) noexcept;
    const elem &operator[](size_t i) const noexcept;
    template<class mylist>
    void init(size_t size, mylist &owner) noexcept;
    void merge(ListStorage &other) noexcept;
    void remove(const size_t &id) noexcept;
    size_t size() const noexcept;
    template<class mylist>
    elem *add(const elem &other, mylist &owner) noexcept;
    typename std::vector<elem *>::iterator begin() noexcept;
    typename std::vector<elem *>::iterator end() noexcept;
    typename std::vector<elem *>::const_iterator begin() const noexcept;
    typename std::vector<elem *>::const_iterator end() const noexcept;
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
    ListStorage(const ListStorage &other, mylist &owner, mylist &storage) noexcept;
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
