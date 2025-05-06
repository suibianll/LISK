#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <bits/hash_bytes.h>
#include "util.h"
#include "hash.h"
#ifndef __CLHT_H__
#define __CLHT_H__
namespace lisk{

#define MAP_INVLD 0
#define MAP_VALID 1
#define MAP_INSRT 2

#define KEY_BUCKT 3
#define ENTRIES_PER_BUCKET KEY_BUCKT
#define CLHT_PERC_EXPANSIONS  1
#define CLHT_MAX_EXPANSIONS   24
#define CLHT_PERC_FULL_DOUBLE 50	   /* % */
#define CLHT_RATIO_DOUBLE     2		  
#define CLHT_OCCUP_AFTER_RES  40
#define CLHT_PERC_FULL_HALVE  5		   /* % */
#define CLHT_RATIO_HALVE      8		  
#define CLHT_MIN_CLHT_SIZE    8

extern uint64_t clht_probe_length;
  /// Round up to next higher power of 2 (return x if it's already a power
  /// of 2) for 32-bit numbers
static inline uint64_t pow2roundup (uint64_t x){
    if (x==0) return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return x+1;
}

union clht_snapshot_t
{
    uint64_t snapshot;
    struct unpakced
    {
    #if KEY_BUCKET == 4
        uint32_t version;
    #elif KEY_BUCKET == 6
        uint16_t version;
    #else
        uint32_t version;
    #endif
        uint8_t map[KEY_BUCKT];
    };
    
};

template<class T, class P>
struct __attribute__ ((aligned (CACHE_LINE_SIZE))) bucket_t
{
    clht_snapshot_t s;
    bucket_t* next;
    T key[KEY_BUCKT];
    P val[KEY_BUCKT];

    inline bool find(const T &k, P &v){
        
        for (size_t i = 0; i < KEY_BUCKT; i++){
            if(k == key[i]){
                v= val[i];
                return true;
            }
        }
        return false;
    }
    inline bool Locate(const T &k, P* &v){
        
        for (size_t i = 0; i < KEY_BUCKT; i++){
            if(k == key[i]){
                v= &val[i];
                return true;
            }
        }
        return false;
    }
};

// std::ofstream out;

template<class T, class P>
class clht
{
private:
    bucket_t<T,P> *buckets;
    uint64_t num_buckets;
    uint64_t hashcode;
    uint32_t num_expands;
    uint32_t num_expands_threshold;
public:
    uint64_t count_;
private:
    char padding[24];
    static inline uint64_t hash(T key, uint64_t h){
        return std::_Hash_bytes(&key, sizeof(T), static_cast<size_t>(0xc70f6907UL)) & h;
    }
    static inline uint64_t hash(T key){
        return std::_Hash_bytes(&key, sizeof(T), static_cast<size_t>(0xc70f6907UL));
    }
public:
    clht(uint64_t num);
    clht(std::pair<T, P> *recs_, uint64_t num);
    clht(T* keys, P* values, uint64_t num);
    ~clht();
    bool Find(const T &key, P &val);
    bool Locate(const T &key, P* &val);
    bool Insert(const T &key, const P &val);
    uint64_t FindInsert(const T &key, const P &val){return 0;}
    bool Remove(const T &key);
    bool Update(const T &key, const P &val){return false;}

    inline bool is_Full(){return false;}//for expanding the capacity
    
    inline bool IsFull(){return false;}//for outer calling

    inline bool IsEmpty(){return count_ == 0;};
    uint64_t Size(){return count_;}
    uint64_t Memory(){ return 0;}
    int resize();
    int do_resize(bool is_increase, int by, std::vector<std::pair<T, P>> &kv_pairs);

    inline Record<T, P> PopFront(){
        Record<T, P> result = {0, 0};
        return result;
    }
    inline Record<T, P> Front(){
        return {0, 0};
    }
    void Dump(_key_t *keys, _value_t* values, uint64_t &start){}
    __always_inline void Prefetch(){
        prefetch(this);
    }
};
template<class T, class P>
clht<T,P>::clht(uint64_t num)
{
    num_buckets = num * 2;
    hashcode = num_buckets - 1;
    buckets =(bucket_t<T,P> *) aligned_alloc(CACHE_LINE_SIZE, num_buckets * sizeof(bucket_t<T,P>));
    memset(buckets, 0, num_buckets * sizeof(bucket_t<T,P>));
    num_expands = 0;
    num_expands_threshold = num_buckets * CLHT_PERC_EXPANSIONS;

    
}
template<class T, class P>
clht<T,P>::clht(std::pair<T, P> *recs_, uint64_t num){
    num_buckets = num * 2;
    hashcode = num_buckets - 1;
    buckets =(bucket_t<T,P> *) aligned_alloc(CACHE_LINE_SIZE, num_buckets * sizeof(bucket_t<T,P>));
    memset(buckets, 0, num_buckets * sizeof(bucket_t<T,P>));
    num_expands = 0;
    num_expands_threshold = num_buckets * CLHT_PERC_EXPANSIONS;

    for (size_t i = 0; i < num; i++){
        Insert(recs_[i].first, recs_[i].second);
    }
}

template<class T, class P>
clht<T,P>::clht(T* keys, P* values, uint64_t num){
    num_buckets = num * 2;
    hashcode = num_buckets - 1;
    buckets =(bucket_t<T,P> *) aligned_alloc(CACHE_LINE_SIZE, num_buckets * sizeof(bucket_t<T,P>));
    memset(buckets, 0, num_buckets * sizeof(bucket_t<T,P>));
    num_expands = 0;
    num_expands_threshold = num_buckets * CLHT_PERC_EXPANSIONS;

    for (size_t i = 0; i < num; i++){
        Insert(keys[i], values[i]);
    }
}

template<class T, class P>
clht<T,P>::~clht()
{
    sizeof(clht<uint64_t, uint64_t>);
    free(buckets);
}
template<class T, class P>
bool clht<T,P>::Find(const T &key, P &val){
#ifndef NDEBUG
    of_search_time++;
#endif
    uint64_t key_hash = hash(key);
    uint64_t bin = key_hash & hashcode;

    bucket_t<T,P> *bucket = &buckets[bin];
    bool found;
    do
    {
        #ifndef NDEBUG
        of_probe_length ++;
        #endif
        bool found = bucket->find(key, val);
        if(found) return true;
        bucket = bucket->next;
    } while (unlikely(bucket != nullptr));
    
    
    return false;
}

template<class T, class P>
bool clht<T,P>::Locate(const T &key, P* &val){
#ifndef NDEBUG
    of_search_time++;
#endif
    uint64_t key_hash = hash(key);

    uint64_t bin = key_hash & hashcode;

    bucket_t<T,P> *bucket = &buckets[bin];
    bool found;
    do
    {
        #ifndef NDEBUG
        of_probe_length ++;
        #endif
        bool found = bucket->Locate(key, val);
        if(found) return true;
        bucket = bucket->next;

    } while (unlikely(bucket != nullptr));
    
    
    return false;
}
template<class T, class P>
bool clht<T,P>::Insert(const T &key, const P &val){
    uint64_t key_hash = hash(key);

    uint64_t bin = key_hash & hashcode;
    bucket_t<T,P> *bucket = &buckets[bin];
    P v;

    if(bucket->find(key, v)){
        return false;
    }

    do
    {
        #ifndef NDEBUG
        clht_probe_length ++;
        #endif
        int empty_pos = -1;

        for (size_t i = 0; i < ENTRIES_PER_BUCKET; i++)
        {
            if(unlikely(bucket->key[i] == key)){//already exists, update value
                bucket->val[i] = val;
                return true;
            }
            else if(bucket->key[i] == 0){
                empty_pos = i;
                break;
            }
        }
        
        if(empty_pos >= 0){
            bucket->key[empty_pos] = key;
            bucket->val[empty_pos] = val;
            return true;
        }
        else{//full
            if(!bucket->next){//no next bucket, need to expand
                bucket_t<T, P>* new_bucket = (bucket_t<T, P> *) aligned_alloc(CACHE_LINE_SIZE, sizeof(bucket_t<T, P>));
                memset(new_bucket, 0 ,sizeof(bucket_t<T, P>));
                new_bucket->key[0] = key;
                new_bucket->val[0] = val;
                bucket->next = new_bucket;
                num_expands++;
                if(num_expands == num_expands_threshold){//expand too many times, need to resize
                    resize();
                }
                return true;
            }
        }
        bucket = bucket->next;
        
    } while (bucket != nullptr);
    
    __builtin_unreachable();
    return false;
}
template<class T, class P>
bool clht<T,P>::Remove(const T &key){
    return false;
}

template<class T, class P>
int clht<T,P>::resize(){
    int num_kv = 0;
    std::vector<std::pair<T, P>> kv_pairs;
    kv_pairs.reserve(100);
    int max_expands = 0;
    for (size_t i = 0; i < num_buckets; i++)//get total number of keys
    {
        int expands_count = -1;
        bucket_t<T, P> *bucket = &buckets[i];
        do
        {
            for (size_t j = 0; j < ENTRIES_PER_BUCKET; j++)
            {
                if(bucket->key[j] != 0){
                    num_kv ++;
                    kv_pairs.emplace_back(std::pair<T,P>(bucket->key[j],bucket->val[j]));
                }
            }
            
            bucket = bucket->next;
            expands_count++;
        } while (bucket);
        if(expands_count > max_expands){
            max_expands = expands_count;
        }
    }

    double fill_ratio = 100.0 * num_kv / (num_buckets * ENTRIES_PER_BUCKET);

    if(fill_ratio > 0 && fill_ratio < CLHT_PERC_FULL_HALVE){//fill ratio < 5%, shrink the hash table by half
        do_resize(false, 0, kv_pairs);
    }
    else if((fill_ratio > CLHT_PERC_FULL_DOUBLE) || (max_expands > CLHT_MAX_EXPANSIONS)){//fill ratio > 50%||too long linked list, expand the hash table
        int inc_by = fill_ratio / CLHT_OCCUP_AFTER_RES;//40% full after resize
        int inc_by_power2 = pow2roundup(inc_by);
        if(inc_by_power2 == 1) inc_by_power2 =2;

        do_resize(true, inc_by_power2, kv_pairs);
    }

    
    return 0;
}

template<class T, class P>
int clht<T,P>::do_resize(bool is_increase, int by, std::vector<std::pair<T, P>> &kv_pairs){
    int new_num_buckets = 0;
    if(is_increase){
        new_num_buckets = num_buckets * by;
    }
    else{
        new_num_buckets = num_buckets / CLHT_RATIO_HALVE;
    }

    bucket_t<T, P> * new_buckets = (bucket_t<T, P> *) aligned_alloc(CACHE_LINE_SIZE, new_num_buckets * sizeof(bucket_t<T, P>));
    memset(new_buckets, 0, new_num_buckets * sizeof(bucket_t<T, P>));

    uint32_t new_num_expands = 0;
    uint32_t new_num_expands_threshold = new_num_buckets * CLHT_PERC_EXPANSIONS;
    uint64_t new_hashcode = new_num_buckets - 1;
    //rehash the keys into the new buckets
    int size = kv_pairs.size();
    for (size_t i = 0; i < size; i++){
        int bin = hash(kv_pairs[i].first, new_hashcode);
        bucket_t<T, P>* bucket = &new_buckets[bin];
        do
        {
            int empty_pos = -1;
            for (size_t j = 0; i < ENTRIES_PER_BUCKET; j++){
                if(bucket->key[j] == 0) empty_pos = j;
            }
            if(empty_pos >= 0){
                bucket->key[empty_pos] = kv_pairs[i].first;
                bucket->val[empty_pos] = kv_pairs[i].second;
                break;
            }
            else if(!bucket->next){//create new bucket
                bucket_t<T, P>* new_bucket = (bucket_t<T, P> *) aligned_alloc(CACHE_LINE_SIZE, sizeof(bucket_t<T, P>));
                memset(new_bucket, 0 ,sizeof(bucket_t<T, P>));
                new_bucket->key[0] = kv_pairs[i].first;
                new_bucket->val[0] = kv_pairs[i].second;
                bucket->next = new_bucket;
                new_num_expands++;
                break;
            }
            bucket = bucket->next;
        } while (true);
    }
    buckets = new_buckets;
    num_expands = new_num_expands;
    num_expands_threshold = new_num_expands_threshold;
    num_buckets = new_num_buckets;
    hashcode = new_hashcode;
        
    if(new_num_expands > new_num_expands_threshold){//expand too many times, need one more resize
        resize();    
    }
    return 0;

}
}//namespace tli

#endif