#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include "util.h"
#include "hash.h"


#define ARRAY_SIZE 7
namespace lisk {

template <typename T, typename P>
class HashNode
{
private:
	uint32_t capacity_;
	uint32_t key_num;
	Record<T, P> *dataslot;
	char buffer[sizeof(Record<T, P>) * ARRAY_SIZE];
private:
	__always_inline size_t key_to_idx(T key){
		return std::_Hash_bytes(&key, sizeof(T), static_cast<size_t>(0xc70f6907UL)) & (capacity_ - 1);
	}
	__always_inline size_t probe_next(int idx) {
		return (idx + 1) & (capacity_ - 1);
	}

	__always_inline int diff(int a, int b){
		return (capacity_ + (a - b)) & (capacity_ - 1);
	}
	uint32_t nextPowerOfTwo(uint32_t x) {
		if (x == 0) return 1; 
		x--; 
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16; 
		return x + 1;
	}

	void emplace(Record<T, P>* data, size_t cap, T key, P val){
		size_t idx = std::_Hash_bytes(&key, sizeof(T), static_cast<size_t>(0xc70f6907UL)) & (cap - 1);
		for (size_t i = idx; ; i = (i + 1) & (cap - 1)){
			if(data[i].key == FREE_FLAG){
				data[i].key = key;
				data[i].val = val;
				return;
			}
		}
	}
	void rehash(int count) {
		count = std::max(count, (int)capacity_ * 2);
		count = nextPowerOfTwo(count);
		Record<T, P>* new_dataslot = new Record<T, P>[count];
		for (size_t i = 0; i < count; i++){
			new_dataslot[i].key = FREE_FLAG;
		}
		
		for (size_t i = 0; i < capacity_; i++){
			if(dataslot[i].key != FREE_FLAG){
				emplace(new_dataslot, count, dataslot[i].key, dataslot[i].val);
			}
		}
		std::swap(dataslot, new_dataslot);
		capacity_ = count;
		if((uint64_t)new_dataslot != (uint64_t) buffer)
			delete[] new_dataslot;
	}

	void reserve(int count) {
		if (count * 2 > capacity_) {
			rehash(count * 2);
		}
	}
public:
	HashNode(T* keys, P* values, int num){
		if(num < ARRAY_SIZE){
			capacity_ = ARRAY_SIZE;
			key_num = num;
			dataslot = reinterpret_cast<Record<T, P> *>(buffer);
			for (size_t i = 0; i < num; i++){
				dataslot[i] = {keys[i], values[i]};
			}
			for (size_t i = num; i < capacity_; i++){
				dataslot[i].key = FREE_FLAG;
			}
		}
		else{
			capacity_ = nextPowerOfTwo(num * 2);
			key_num = num;
			dataslot = new Record<T, P>[capacity_];
			for (size_t i = 0; i < capacity_; i++){
				dataslot[i].key = FREE_FLAG;
			}
			for (size_t i = 0; i < num; i++){
				emplace(dataslot, capacity_, keys[i], values[i]);
			}
		}
	}
	HashNode(int num){
		capacity_ = ARRAY_SIZE;
		key_num = 0;
		dataslot = reinterpret_cast<Record<T, P> *>(buffer);
		for (size_t i = 0; i < capacity_; i++){
			dataslot[i].key = FREE_FLAG;
		}
	}
	~HashNode(){delete[] dataslot;}

	inline bool is_Full(){return key_num == capacity_;}//for expanding the capacity
    
    inline bool IsFull(){return key_num == 12 << 6;}//for outer calling

    inline bool IsEmpty(){return key_num == 0;};

    uint64_t Size(){
        return key_num;
    }

    uint64_t Memory(){
        return sizeof(*this) + sizeof(Record<T, P>) * capacity_;
    }
	bool Find(const T &key, P &val){
		#ifndef NDEBUG
            of_search_time++;
        #endif
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				#ifndef NDEBUG
                    of_probe_length++;
                #endif
				if(dataslot[i].key == key){
					val = dataslot[i].val;
					return true;
				}
			}
			return false;
		}
		for (size_t i = key_to_idx(key);; i = probe_next(i)){
			#ifndef NDEBUG
                of_probe_length++;
            #endif
			if(dataslot[i].key == key){
				val = dataslot[i].val;
				return true;
			}
			if(dataslot[i].key == FREE_FLAG) return false;
		}
		return false;
	}
	bool Locate(const T &key, P* &val){
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				#ifndef NDEBUG
                    of_probe_length++;
                #endif
				if(dataslot[i].key == key){
					val = &dataslot[i].val;
					return true;
				}
			}
			return false;
		}
		for (size_t i = key_to_idx(key);; i = probe_next(i)){
			#ifndef NDEBUG
                of_probe_length++;
            #endif
			if(dataslot[i].key == key){
				val = &dataslot[i].val;
				return true;
			}
			if(dataslot[i].key == FREE_FLAG) return false;
		}
		return false;
	}

	uint64_t FindInsert(const T &key, P &val){
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				#ifndef NDEBUG
                    of_probe_length++;
                #endif
				if(dataslot[i].key == key){
					val = dataslot[i].val;
					return true;
				}
			}
			return false;
		}
		for (size_t i = key_to_idx(key);; i = probe_next(i)){
			if(dataslot[i].key == key){
				val = dataslot[i].val;
				return uint64_t(&dataslot[i].val);
			}
			if(dataslot[i].key == FREE_FLAG){
				dataslot[i].key = key;
				dataslot[i].val = val;
				return 0;
			}
		}
	}
	bool Insert(const T &key, const P &val){
		if(key_num < ARRAY_SIZE){
			size_t i = 0;
			for (; i < key_num; i++){
				if(dataslot[i].key > key) break;
				else if(dataslot[i].key == key){
					return false;
				}
			}
			memmove(&dataslot[i + 1], &dataslot[i], sizeof(Record<T,P>) * (key_num - i));
			dataslot[i] = Record<T,P>(key, val);
			key_num ++;
			if(key_num == ARRAY_SIZE){
				rehash(capacity_ * 2);
			}
			return true;
		}
		reserve(key_num + 1);
		for (size_t i = key_to_idx(key); ; i = probe_next(i)){
			if(dataslot[i].key == FREE_FLAG){
				dataslot[i].key = key;
				dataslot[i].val = val;
				key_num ++;
				return true;
			}
			else if(dataslot[i].key == key){
				return false;
			}
		}
		return false;
	}
	bool Update(const T &key, const P &val){
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				if(dataslot[i].key == key){
					dataslot[i].val = val;
					return true;
				}
			}
			return false;
		}
		for (size_t i = key_to_idx(key);; i = probe_next(i)){
			if(dataslot[i].key == key){
				dataslot[i].val = val;
				return true;
			}
			if(dataslot[i].key == FREE_FLAG) return false;
		}
		return false;
	}
	bool Remove(const T &key){
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				if(dataslot[i].key == key){
					dataslot[i].key = FREE_FLAG;
					return true;
				}
			}
			return false;
		}
		size_t bucket = key_to_idx(key);
		for(;;bucket = probe_next(bucket)){
			if(dataslot[bucket].key == key) break;
			if(dataslot[bucket].key == FREE_FLAG) return false;
		}
		for (size_t i = probe_next(bucket);; i = probe_next(i)){
			if(dataslot[i].key == FREE_FLAG){
				dataslot[bucket].key = FREE_FLAG;
				key_num --;
				return true;
			}
			size_t ideal = key_to_idx(dataslot[i].key);
			if(diff(bucket, ideal) < diff(i, ideal)){
				dataslot[bucket] = dataslot[i];
				bucket = i;
			}
		}
		
		return true;
	}
	int Scan(T &startKey, int len, Record<T,P> *result){return 0;}

	void Prefetch(){prefetch(this);}//prefetch(this);

	void Dump(_key_t *keys, _value_t* values, uint64_t &start){
		if(key_num < ARRAY_SIZE){
			for (size_t i = 0; i < key_num; i++){
				keys[start] = dataslot[i].key;
				values[start] = dataslot[i].val;
				start ++;
			}
		}
		else{
			int count = 0;
			std::sort(dataslot, dataslot + capacity_);
			for (size_t i = 0; i < key_num; i++){
				keys[start + count] = dataslot[i].key;
				values[start + count] = dataslot[i].val;
				count ++;
			}
			start += count;
		}
		
		
	}

	inline Record<T,P> PopFront(){return Record<T,P>();}

	inline Record<T,P> Front(){return Record<T,P>();}

};


} // namespace lisk