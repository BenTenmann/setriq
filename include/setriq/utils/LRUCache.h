//
// Taken from https://github.com/lamerman/cpp-lru-cache/blob/master/include/lrucache.hpp
//

#ifndef METRICS_LRUCACHE_H
#define METRICS_LRUCACHE_H

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>

template<typename key_t, typename value_t>
class LRUCache {
    typedef typename std::pair<key_t, value_t> keyValuePairT;
    typedef typename std::list<keyValuePairT>::iterator listIteratorT;

private:
    std::list<keyValuePairT> cacheItemsList;
    std::unordered_map<key_t, listIteratorT> cacheItemsMap;
    size_t maxSize;

public:
    LRUCache() : maxSize {0} {};
    LRUCache(size_t size) : maxSize {size} {};

    void put(const key_t& key, const value_t& value) {
        auto it = this->cacheItemsMap.find(key);
        this->cacheItemsList.push_front(keyValuePairT(key, value));
        if (it != this->cacheItemsMap.end()) {
            this->cacheItemsList.erase(it->second);
            this->cacheItemsMap.erase(it);
        }
        this->cacheItemsMap[key] = this->cacheItemsList.begin();

        if (this->cacheItemsMap.size() > maxSize) {
            auto last = this->cacheItemsList.end();
            last--;
            this->cacheItemsMap.erase(last->first);
            this->cacheItemsList.pop_back();
        }
    }

    const value_t& get(const key_t& key) {
        auto it = this->cacheItemsMap.find(key);
        if (it == this->cacheItemsMap.end()) {
            throw std::range_error("There is no such key in cache");
        } else {
            this->cacheItemsList.splice(this->cacheItemsList.begin(), this->cacheItemsList, it->second);
            return it->second->second;
        }
    }

    bool exists(const key_t& key) const {
        return this->cacheItemsMap.find(key) != this->cacheItemsMap.end();
    }

    size_t size() const {
        return this->cacheItemsMap.size();
    }
};

#endif //METRICS_LRUCACHE_H
