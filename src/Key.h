/*** 
 * @Author: Chu Zhaole
 * @Date: 2024-01-17 02:25:10
 * @LastEditTime: 2025-05-06 08:38:23
 * @LastEditors: Chu Zhaole
 * @Description: 
 * @FilePath: /LISK/src/Key.h
 * @Copyright (c) Chu Zhaole.
 */

#ifndef __LISK_KEY_H__
#define __LISK_KEY_H__

#include <stdint.h>
#include <cstring>
#include <memory>
#include <assert.h>
// #include "util.h"

// #define ENCODE
namespace lisk{
//key type = 0(uint64), 1(double)
#define KEYTYPE 0
using _key_t = uint64_t;
using _value_t = uint64_t;
using KeyLen = uint32_t;
static const uint64_t SLICE = 8;
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
class Key {
    static constexpr uint32_t stackLen = 40;
public:
    uint32_t len = 0;

    _value_t val;

    uint8_t *data;

    uint8_t stackKey[stackLen];


    Key() {len = 0;}

    ~Key();

    Key(const Key &key);

    Key(Key &&key);

    Key(uint64_t k) { setInt(k); }

    void set(const char bytes[], const std::size_t length);

    void set(const char bytes[], const std::size_t length, const _value_t v);

    __always_inline void setInt(uint64_t k) { data = stackKey; len = 8; *reinterpret_cast<uint64_t*>(stackKey) = k; }
    __always_inline void setInt(uint64_t k, uint32_t depth){
        *reinterpret_cast<uint64_t*>(data + depth * 8) = (k);
    }

    __always_inline void setVal(_value_t v) {val = v;}

    __always_inline _key_t getSlice(uint32_t depth) const {
        _key_t  result = 0;
#if KEYTYPE == 0
        if(likely(depth + SLICE <= len)){
            result = *reinterpret_cast<_key_t*>(data + depth);
            return __bswap_64(result);//reverse byte
        }
        else{
            for (size_t i = depth; i < len; i++){
                result <<= 8;
                result += data[i];
            }
            result <<= (SLICE - (len - depth))*8; 
            return result;   
        }
#elif KEYTYPE == 1
        double weight = 1.0;
        for(int i = depth; i < depth + SLICE && i < len; i++){
            result += (double)(data[i]) * weight;
            weight /= 256;
        }
        return result;
#endif
        
    }

    void operator=(const char key[]);
    void operator=(const Key &k);

    inline int keycompare(const Key &other) const{
        int len = this->len < other.len ? this->len: other.len;
        int cmp = memcmp(data, other.data, len);
        if(cmp == 0){
            if(this->len < other.len)
                return -1;
            else if(this->len > other.len)
                return 1;
            else
                return 0;
        }
        return cmp;
    }
    inline bool operator==(const Key &other) const { 
        if(len != other.len)
            return false;
        return memcmp(data, other.data, len) == 0;
    }
    inline bool operator!=(const Key &other) const { return !(*this == other);}

    inline bool operator<(const Key &other) const { return keycompare(other) < 0;}
    inline bool operator>(const Key &other) const { return keycompare(other) > 0;}
    inline bool operator<=(const Key &other) const { return !(*this > other);}
    inline bool operator>=(const Key &other) const {return !(*this < other);}

    
    inline int keycompare(const Key &other, uint16_t depth) const{
        int len = this->len < other.len ? this->len: other.len;
        int cmp = memcmp(data + depth, other.data + depth, len - depth);
        if(cmp == 0){
            if(this->len < other.len)
                return -1;
            else if(this->len > other.len)
                return 1;
            else
                return 0;
        }
        return cmp;
    }
    inline bool verify(Key &other, uint16_t depth) const{
        if(len != other.len) return false;
        if(depth >= len) return true;
        return memcmp(data + depth, other.data + depth, len - depth) == 0;
    }
    uint8_t &operator[](std::size_t i);

    const uint8_t &operator[](std::size_t i) const;

    KeyLen getKeyLen() const;

    void setKeyLen(KeyLen len);

};


__always_inline uint8_t &Key::operator[](std::size_t i) {
    // assert(i < len);
    return data[i];
}

__always_inline const uint8_t &Key::operator[](std::size_t i) const {
    // assert(i < len);
    return data[i];
}

__always_inline KeyLen Key::getKeyLen() const { return len; }

inline Key::~Key() {

}

inline Key::Key(const Key &key) {
    len = key.len;
    val = key.val;
    if (len >= stackLen) {
        data = key.data;
    } else {
        memcpy(stackKey, key.stackKey, key.len);
        data = stackKey;
    }
    data[len] = 0;
}

inline Key::Key(Key &&key) {
    len = key.len;
    val = key.val;
    if (len >= stackLen) {
        data = key.data;
        key.data = nullptr;
    } else {
        memcpy(stackKey, key.stackKey, key.len);
        data = stackKey;
    }
    data[len] = 0;
}

inline void Key::set(const char bytes[], const std::size_t length) {
    if (len >= stackLen) {
        delete[] data;
    }
    if (length < stackLen) {
        memcpy(stackKey, bytes, length);
        data = stackKey;
    } else {
        data = new uint8_t[length+1];
        memcpy(data, bytes, length);
    }
    len = length;
    data[len] = 0;
}

inline void Key::set(const char bytes[], const std::size_t length, const _value_t v) {
    if (len >= stackLen) {
        delete[] data;
    }
    if (length < stackLen) {
        memcpy(stackKey, bytes, length);
        data = stackKey;
    } else {
        data = new uint8_t[length+1];
        memcpy(data, bytes, length);
    }
    val = v;
    len = length;
    data[len] = 0;
}

inline void Key::operator=(const char key[]) {
    if (len >= stackLen) {
        delete[] data;
    }
    len = strlen(key);
    if (len < stackLen) {
        memcpy(stackKey, key, len);
        data = stackKey;
    } else {
        data = new uint8_t[len+1];
        memcpy(data, key, len);
    }
    data[len] = 0;
}

inline void Key::operator=(const Key &k) {
    if (len >= stackLen) {
        delete[] data;
    }
    len = k.getKeyLen();
    if (len < stackLen) {
        memcpy(stackKey, k.data, len);
        data = stackKey;
    } else {
        data = new uint8_t[len+1];
        memcpy(data, k.data, len);
    }
    data[len] = 0;
}

inline void Key::setKeyLen(KeyLen newLen) {
    if (len == newLen) return;
    if (len > stackLen) {
        delete[] data;
    }
    len = newLen;
    if (len > stackLen) {
        data = new uint8_t[len];
    } else {
        data = stackKey;
    }
}

inline int commonPrefixLen(Key &k1, Key &k2){
    int i = 0;
    for(; i < k1.len && i < k2.len && k1[i] == k2[i]; i++){}
    
    return i;
}

__always_inline uint16_t hashKey(const Key &key){
    uint16_t ret = key.len;
    uint16_t c1 = key[ret / 2];
    uint16_t c2 = key[2 * ret / 3];
    uint16_t c3 = key[4 * ret / 5];
    return ret ^ c1 ^ c2 ^ c3;
}

__always_inline Key* setHashKey(const Key *key){
    uint64_t hash = hashKey(*key);
    uint64_t ptr = (uint64_t) (void *)key;
    ptr = ptr | (hash << 48);
    return (Key*) (void *)ptr;
}

__always_inline uint16_t getHash(const Key *key){
    uint64_t v = (uint64_t)(void *)key;
    return (v >> 48) & 0xffff;
}
}
#endif // LIVAK_KEY_H