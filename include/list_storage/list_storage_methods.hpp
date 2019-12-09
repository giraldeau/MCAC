#ifndef INCLUDE_LIST_STORAGE_LIST_STORAGE_METHODS_HPP
#define INCLUDE_LIST_STORAGE_LIST_STORAGE_METHODS_HPP 1
#include "list_storage/list_storage.hpp"
#include "physical_model/physical_model.hpp"
#include <array>
#include <vector>


namespace mcac {
template<int N, class elem>
[[gnu::pure]] size_t ListStorage<N, elem>::size() const noexcept {
    return list.size();
}
template<int N, class elem>
[[gnu::pure]] elem &ListStorage<N, elem>::operator[](size_t i) noexcept {
    return *list[i];
}
template<int N, class elem>
[[gnu::pure]] const elem &ListStorage<N, elem>::operator[](size_t i) const noexcept {
    return *list[i];
}
template<int N, class elem>
void ListStorage<N, elem>::destroy() noexcept {
    if (!external_storage) {
        if (storage) {
            if (size() > 0) {
                for (size_t n = size(); n-- > 0;) {
                    delete list[n];
                }
            }
            delete storage;
        }
    } else {
        storage = nullptr;
    }
}
//template<int N, class elem>
//void swap(ListStorage<N, elem> &first, ListStorage<N, elem> &second) noexcept {
//    using std::swap;
//    swap(first.list, second.list);
//    swap(first.storage, second.storage);
//    swap(first.external_storage, second.external_storage);
//}
template<int N, class elem>
void ListStorage<N, elem>::merge(ListStorage<N, elem> &other) noexcept {
    if (!external_storage) {
        for (size_t i = 0; i < N; i++) {
            (*storage)[i].insert((*storage)[i].end(), (*other.storage)[i].begin(), (*other.storage)[i].end());
        }
        for (auto &data : other.list) {
            list.push_back(data);
            data = nullptr;
        }
    } else {
        list.insert(list.end(), other.list.begin(), other.list.end());
    }
}
template<int N, class elem>
void ListStorage<N, elem>::remove(const size_t &id) noexcept {
    delete list[id];
    list.erase(list.begin() + long(id));
    for (size_t i = 0; i < N; i++) {
        (*storage)[i].erase((*storage)[i].begin() + long(id));
    }
    for (size_t i = id; i < list.size(); i++) {
        list[i]->decrease_label();
    }
}
template<int N, class elem>
typename std::vector<elem *>::iterator ListStorage<N, elem>::begin() noexcept {
    return list.begin();
}
template<int N, class elem>
typename std::vector<elem *>::iterator ListStorage<N, elem>::end() noexcept {
    return list.end();
}
template<int N, class elem>
typename std::vector<elem *>::const_iterator ListStorage<N, elem>::begin() const noexcept {
    return list.begin();
}
template<int N, class elem>
typename std::vector<elem *>::const_iterator ListStorage<N, elem>::end() const noexcept {
    return list.end();
}
template<int N, class elem>
template<class mylist>
elem *ListStorage<N, elem>::add(const elem &other, mylist &owner) noexcept {
    return new elem(other, &owner);
}
/** Default constructor in local storage */
//template<int N, class elem>
//template<class mylist>
//ListStorage<N, elem>::ListStorage(mylist &owner, const PhysicalModel &the_physical_model) noexcept:
//    list(the_physical_model.N),
//    storage(new std::array<std::vector<double>, N>),
//    //    storage(nullptr),
//    external_storage(nullptr) {
//    for (std::vector<double> &data : (*storage)) {
//        data.assign(the_physical_model.N, 0.);
//    }
//    for (size_t i = 0; i < the_physical_model.N; i++) {
//        list[i] = new elem(owner, i, the_physical_model);
//    }
//}
///** Default constructor in local storage */
template<int N, class elem>
ListStorage<N, elem>::ListStorage() noexcept:
    list(),
    storage(nullptr),
    external_storage(nullptr) {
}
template<int N, class elem>
template<class mylist>
void ListStorage<N, elem>::init(size_t size, mylist &owner) noexcept {
    destroy();
    storage = new std::array<std::vector<double>, N>;
    external_storage = nullptr;

    //preallocation
    list.reserve(size);

    //initialize data
    for (std::vector<double> &data : (*storage)) {
        data.assign(size, 0.);
    }
    const size_t _list_size = size;
    for (size_t i = 0; i < _list_size; i++) {
        list.push_back(new elem(&owner, i));
    }
}
/** Constructor with external storage */
template<int N, class elem>
ListStorage<N, elem>::ListStorage(ListStorage<N, elem> &parent, const std::vector<size_t> &index) noexcept:
    list(),
    storage(parent.storage),
    external_storage(&parent) {
    list.assign(index.size(), nullptr);
    const size_t _list_size = size();
//    #pragma omp for simd
        for (size_t i = 0; i < _list_size; i++) {
            list[i] = external_storage->list[index[i]];
        }
}
/** Destructor */
template<int N, class elem>
ListStorage<N, elem>::~ListStorage() noexcept {
    destroy();
}
/** Copy constructor */
template<int N, class elem>
template<class mylist>
ListStorage<N, elem>::ListStorage(const ListStorage<N, elem> &other, mylist &owner, mylist &ext_storage) noexcept:
    list(),
    storage(ext_storage.storage),
    external_storage(&ext_storage) {
    const size_t _start = (*storage)[0].size();
    for (size_t i = 0; i < N; i++) {
        (*storage)[i].reserve(other.size() + _start);
    }
    list.reserve(other.size() + _start);
    const size_t _list_size = other.size();
    for (size_t size = 0; size < _list_size; size++) {
        list.push_back(new elem(*other.list[size], &ext_storage, size + _start));
    }
}
/** Move constructor */
template<int N, class elem>
ListStorage<N, elem>::ListStorage(ListStorage &&other) noexcept :
    list(other.list),
    storage(other.storage),
    external_storage(other.external_storage) {
    storage = nullptr;
    external_storage = nullptr;
}
/** Move assignment operator */
template<int N, class elem>
ListStorage<N, elem> &ListStorage<N, elem>::operator=(ListStorage<N, elem> &&other) noexcept {
    destroy();
    storage = other.storage;
    external_storage = other.external_storage;
    list = std::move(other.list);
    other.storage = nullptr;
    other.external_storage = nullptr;
    return *this;
}
/** Copy constructor */
//template<int N, class elem>
//template<class mylist>
//ListStorage<N, elem>::ListStorage(const ListStorage<N, elem> &other, const mylist &owner) noexcept:
//    list(),
//    storage(new std::array<std::vector<double>, N>),
//    external_storage(nullptr) {
//    for (size_t i = 0; i < N; i++) {
//        (*storage)[i].reserve(other.size());
//    }
//    const size_t listSize = other.size();
//    for (size_t _size = 0; _size < listSize; _size++) {
//        for (size_t i = 0; i < N; i++) {
//            (*storage)[i][_size] = (*other.storage)[i][other.list[_size]->indexInStorage];
//        }
//        list.push_back(new elem(owner, _size));
//    }
//}
///** Copy assignment operator */
//template<int N, class elem>
//ListStorage<N, elem> &ListStorage<N, elem>::operator=(const ListStorage<N, elem> &other) noexcept {
//    ListStorage<N, elem> tmp{other};         // re-use copy-constructor
//    *this = std::move(tmp);                  // re-use move-assignment
//    return *this;
//}
}  // namespace mcac
#endif //INCLUDE_LIST_STORAGE_LIST_STORAGE_METHODS_HPP
