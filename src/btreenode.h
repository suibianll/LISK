#ifndef __LISK_BTREENODE_H__
#define __LISK_BTREENODE_H__
#include <cassert>
#include <cstring>
#include <atomic>
#include <immintrin.h>
#include <sched.h>
#include "node.h"
#include "util.h"

namespace lisk {
static const uint64_t pageSize=256;

class BTreeLeafNode: public Node{
public:
static const uint64_t maxEntries=(pageSize-sizeof(Node)-sizeof(BTreeLeafNode *))/(sizeof(Record<_key_t, _value_t>));

BTreeLeafNode* sibling_;
Record<_key_t, _value_t> dataslot[maxEntries];

public:
	BTreeLeafNode(uint16_t d = 0){
		depth = d;
		key_num = 0;
		node_type = BLEAF;
		prefix_len = 0;
		sibling_ = nullptr;
	}
	BTreeLeafNode(_key_t *keys, _value_t *values, uint64_t num, uint16_t d = 0){
		depth = d;
		key_num = num;
		node_type = BLEAF;
		prefix_len = 0;
		sibling_ = nullptr;
		for (size_t i = 0; i < num; i++){
			dataslot[i].key = keys[i];
			dataslot[i].val = values[i];
		}
		return;
	}
	void bulk_load(_key_t *keys, _value_t *values, uint64_t num, uint16_t d = 0){
		depth = d;
		key_num = num;
		for (size_t i = 0; i < num; i++){
			dataslot[i].key = keys[i];
			dataslot[i].val = values[i];
		}
		return;
	}
	bool isFull() { return key_num==maxEntries; };

    unsigned lowerBound(_key_t k) {
        unsigned lower = 0;
        unsigned upper = key_num;
        do {
			unsigned mid = ((upper - lower) / 2) + lower;
			if (k < dataslot[mid].key) {
				upper=mid;
			} else if (k > dataslot[mid].key) {
				lower=mid+1;
			} else {
				return mid;
			}
        } while (lower < upper);
        return lower;
    }

    void insert(_key_t k, _value_t p) {
      if (key_num) {
        unsigned pos=lowerBound(k);
        if ((pos<key_num) && (dataslot[pos].key == k)) {
			// Upsert
			dataslot[pos].val = p;
			return;
        }
		memmove(dataslot + pos + 1, dataslot + pos, sizeof(Record<_key_t, _value_t>) * (key_num - pos));
        dataslot[pos].key = k;
        dataslot[pos].val = p;
      } else {
        dataslot[0].key = k;
        dataslot[0].val = p;
      }
      key_num++;
    }

	void append(_key_t k, _value_t p){
        dataslot[key_num].key = k;
        dataslot[key_num].val = p;
		key_num ++;
	}

    BTreeLeafNode* split(_key_t& sep) {
        BTreeLeafNode* newLeaf = new BTreeLeafNode();
        newLeaf->key_num = key_num - (key_num / 2);
        key_num = key_num - newLeaf->key_num;
		memcpy(newLeaf->dataslot, dataslot + key_num, sizeof(Record<_key_t, _value_t>) * newLeaf->key_num);
        newLeaf->sibling_ = sibling_;
        sibling_ = newLeaf;
        sep = dataslot[key_num - 1].key;
        return newLeaf;
    }

	uint64_t Memory() const{
		return pageSize;
	}
};

class BTreeInnerNode: public Node{
public:
	static const uint64_t maxEntries=(pageSize-sizeof(Node))/(sizeof(Record<_key_t, Node*>));

	Record<_key_t, Node*> datatslot[maxEntries];

public:

    BTreeInnerNode(uint16_t d) {
        depth = d;
		key_num = 0;
		node_type = BINNER;
		prefix_len = 0;
    }

    bool isFull() { return key_num==(maxEntries-1); };

    unsigned lowerBound(_key_t k) {
        unsigned lower = 0;
        unsigned upper = key_num;
        do {
			unsigned mid = ((upper - lower) / 2) + lower;
			if (k < datatslot[mid].key) {
				upper = mid;
			} else if (k > datatslot[mid].key) {
				lower = mid + 1;
			} else {
				return mid;
			}
        } while (lower < upper);
        return lower;
    }

    BTreeInnerNode* split(_key_t& sep) {
        BTreeInnerNode* newInner = new BTreeInnerNode(depth);
        newInner->key_num = key_num - (key_num / 2);
        key_num = key_num-newInner->key_num - 1;
        sep = datatslot[key_num].key;
		memcpy(newInner->datatslot, datatslot + key_num + 1, sizeof(Record<_key_t, Node*>) * (newInner->key_num + 1));
        return newInner;
    }

    void insert(_key_t k, Node* child) {
        assert(key_num<maxEntries-1);
        unsigned pos=lowerBound(k);
		memmove(datatslot + pos + 1, datatslot + pos, sizeof(Record<_key_t, Node*>) * (key_num - pos + 1));
		datatslot[pos].key = k;
		datatslot[pos].val = child;
        std::swap(datatslot[pos].val, datatslot[pos + 1].val);
        key_num++;
    }
	inline void append(_key_t k, Node* child) {
        assert(key_num<maxEntries-1);
		datatslot[key_num].key = k;
		datatslot[key_num].val = child;
        key_num++;
    }

	inline void append(Node* child) {
        assert(key_num<maxEntries-1);
		datatslot[key_num].val = child;
    }
	uint64_t Memory() const{
		return pageSize;
	}
};

extern Node* BTree(_key_t *keys, _value_t *values,  uint64_t num, uint16_t depth);

extern bool BTreeFind(Node* root, Key &key, _value_t &v, int depth);
extern Node* BTreeInsert(Node* root, Key &key, _value_t v, int depth);
extern bool BTreeUpdate(Node* root, Key &key, _value_t v, int depth);
extern bool BTreeRemove(Node* root, Key &key, int depth);
extern uint64_t BTreeScan(Key &key, int range, Record<_key_t, _value_t>* output, int depth);

}//namespace lisk
#endif