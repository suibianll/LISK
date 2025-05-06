#ifndef __LISK_BTREE_H__
#define __LISK_BTREE_H__
#include <cassert>
#include <cstring>
#include <atomic>
#include <immintrin.h>
#include <sched.h>
#include "util.h"

namespace lisk {
extern uint64_t of_probe_length;
namespace btree{
enum class PageType : uint16_t { BTreeInner=1, BTreeLeaf=2 };

static const uint64_t pageSize=256;

struct NodeBase{
  PageType type;
  uint16_t count;
}__attribute((aligned(8)));

struct BTreeLeafBase : public NodeBase {
   static const PageType typeMarker=PageType::BTreeLeaf;
};

template<class Key,class Payload>
struct BTreeLeaf : public BTreeLeafBase {
    struct Entry {
        Key k;
        Payload p;
    };

    static const uint64_t maxEntries=(pageSize-sizeof(NodeBase)-sizeof(BTreeLeaf<Key, Payload> *))/(sizeof(Key)+sizeof(Payload));
    
    // Singly linked list pointer to my sibling
    BTreeLeaf<Key, Payload> *next_leaf;

    Key keys[maxEntries];
    Payload payloads[maxEntries];

    BTreeLeaf() {
        count=0;
        type=typeMarker;
        next_leaf=nullptr;
    }

	void bulkload(std::pair<Key, Payload>* recs_, uint64_t num){
		for (size_t i = 0; i < num; i++){
			keys[i] = recs_[i].first;
			payloads[i] = recs_[i].second;
		}
		return;
	}

    bool isFull() { return count==maxEntries; };

    unsigned lowerBound(Key k) {
        unsigned lower=0;
        unsigned upper=count;
        do {
			#ifndef NDEBUG
			of_probe_length++;
			#endif
			unsigned mid=((upper-lower)/2)+lower;
			if (k<keys[mid]) {
				upper=mid;
			} else if (k>keys[mid]) {
				lower=mid+1;
			} else {
				return mid;
			}
        } while (lower<upper);
        return lower;
    }

    unsigned lowerBoundBF(Key k) {
        auto base=keys;
        unsigned n=count;
        while (n>1) {
			const unsigned half=n/2;
			base=(base[half]<k)?(base+half):base;
			n-=half;
        }
        return (*base<k)+base-keys;
    }

    void insert(Key k,Payload p) {
      assert(count<maxEntries);
      if (count) {
        unsigned pos=lowerBound(k);
        if ((pos<count) && (keys[pos]==k)) {
			// Upsert
			payloads[pos] = p;
			return;
        }
        memmove(keys+pos+1,keys+pos,sizeof(Key)*(count-pos));
        memmove(payloads+pos+1,payloads+pos,sizeof(Payload)*(count-pos));
        keys[pos]=k;
        payloads[pos]=p;
      } else {
        keys[0]=k;
        payloads[0]=p;
      }
      count++;
    }

	void append(Key k, Payload p){
		assert(count<maxEntries);
		keys[count] = k;
		payloads[count] = p;
		count++;
	}

    BTreeLeaf* split(Key& sep) {
        BTreeLeaf* newLeaf = new BTreeLeaf();
        newLeaf->count = count-(count/2);
        count = count-newLeaf->count;
        memcpy(newLeaf->keys, keys+count, sizeof(Key)*newLeaf->count);
        memcpy(newLeaf->payloads, payloads+count, sizeof(Payload)*newLeaf->count);
        newLeaf->next_leaf = next_leaf;
        next_leaf = newLeaf;
        sep = keys[count-1];
        return newLeaf;
    }
  };

  struct BTreeInnerBase : public NodeBase {
    static const PageType typeMarker=PageType::BTreeInner;
  };

  template<class Key>
  struct BTreeInner : public BTreeInnerBase {
    static const uint64_t maxEntries=(pageSize-sizeof(NodeBase))/(sizeof(Key)+sizeof(NodeBase*));
    NodeBase* children[maxEntries];
    Key keys[maxEntries];

    BTreeInner() {
        count=0;
        type=typeMarker;
    }

    bool isFull() { return count==(maxEntries-1); };

    unsigned lowerBoundBF(Key k) {
        auto base=keys;
        unsigned n=count;
        while (n>1) {
			const unsigned half=n/2;
			base=(base[half]<k)?(base+half):base;
			n-=half;
        }
        return (*base<k)+base-keys;
    }

    unsigned lowerBound(Key k) {
        unsigned lower=0;
        unsigned upper=count;
        do {
			#ifndef NDEBUG
			of_probe_length++;
			#endif
			unsigned mid=((upper-lower)/2)+lower;
			if (k<keys[mid]) {
					upper=mid;
			} else if (k>keys[mid]) {
					lower=mid+1;
			} else {
					return mid;
			}
        } while (lower<upper);
        return lower;
    }

    BTreeInner* split(Key& sep) {
        BTreeInner* newInner=new BTreeInner();
        newInner->count=count-(count/2);
        count=count-newInner->count-1;
        sep=keys[count];
        memcpy(newInner->keys,keys+count+1,sizeof(Key)*(newInner->count+1));
        memcpy(newInner->children,children+count+1,sizeof(NodeBase*)*(newInner->count+1));
        return newInner;
    }

    void insert(Key k,NodeBase* child) {
        assert(count<maxEntries-1);
        unsigned pos=lowerBound(k);
        memmove(keys+pos+1,keys+pos,sizeof(Key)*(count-pos+1));
        memmove(children+pos+1,children+pos,sizeof(NodeBase*)*(count-pos+1));
        keys[pos]=k;
        children[pos]=child;
        std::swap(children[pos],children[pos+1]);
        count++;
    }
	inline void append(Key k,NodeBase* child) {
        assert(count<maxEntries-1);
        keys[count]=k;
        children[count]=child;
        count++;
    }

	inline void append(NodeBase* child) {
        assert(count<maxEntries-1);
        children[count]=child;
    }

  };


  template<class Key,class Value>
  struct BTree {
    NodeBase* root;
	BTreeLeaf<Key, Value>* first_leaf;
	uint64_t count_;
    BTree(uint64_t cap = 0) {
        root = new BTreeLeaf<Key,Value>();
		first_leaf = static_cast<BTreeLeaf<Key,Value>*>(root);
    }

	BTree(Key *keys, Value *values, int num){
		Bulkload(keys, values, num);
		// root = new BTreeLeaf<Key,Value>();
		first_leaf = static_cast<BTreeLeaf<Key,Value>*>(root);
		// for (size_t i = 0; i < num; i++)
		// {
		// 	// if(i == 32387)
        //     // 	int ss = 0;
		// 	Insert(recs_[i].first, recs_[i].second);
		// }
		
	}

	BTree(std::pair<Key, Value> * recs_, int num){
		Bulkload(recs_, num);
		// root = new BTreeLeaf<Key,Value>();
		first_leaf = static_cast<BTreeLeaf<Key,Value>*>(root);
		// for (size_t i = 0; i < num; i++)
		// {
		// 	// if(i == 32387)
        //     // 	int ss = 0;
		// 	Insert(recs_[i].first, recs_[i].second);
		// }
		
	}

	//实现bulkload，参数变成(Key *keys, Value *values, int num)
	void Bulkload(Key *keys, Value *values, int num) {
		uint64_t num_leaves = (num + BTreeLeaf<Key, Value>::maxEntries - 1) / BTreeLeaf<Key, Value>::maxEntries;
		uint64_t num_items = num;
	
		std::pair<Key, uint64_t> *nextlevel = new std::pair<Key, uint64_t>[num_leaves];
	
		BTreeLeaf<Key, Value>* last_leaf = nullptr;
		uint64_t cur = 0;
		for (size_t i = 0; i < num_leaves; ++i) {
			BTreeLeaf<Key, Value>* leaf = new BTreeLeaf<Key, Value>();
			uint64_t bulk_num = num_items / (num_leaves - i);
			for (size_t j = 0; j < bulk_num; ++j) {
				leaf->append(keys[cur], values[cur]);
				++cur;
			}
			nextlevel[i].first = keys[cur - 1];
			nextlevel[i].second = (uint64_t)leaf;
			if (last_leaf != nullptr) {
				last_leaf->next_leaf = leaf;
			} else {
				first_leaf = leaf;
			}
			last_leaf = leaf;
	
			num_items -= bulk_num;
		}
	
		if (first_leaf == last_leaf) {
			root = first_leaf;
			delete[] nextlevel;
			return;
		}
	
		// Build inner nodes
		uint64_t num_parents = num_leaves;
		while (num_parents != 1) {
			uint64_t num_children = num_parents;
			num_parents = (num_children + BTreeInner<Key>::maxEntries - 1) / (BTreeInner<Key>::maxEntries - 1);
			uint64_t inner_index = 0;
			for (size_t i = 0; i < num_parents; i++) {
				uint64_t bulk_num = num_children / (num_parents - i) - 1;
				BTreeInner<Key>* new_inner = new BTreeInner<Key>();
				for (size_t j = 0; j < bulk_num; j++) {
					new_inner->append(nextlevel[inner_index].first, (NodeBase *)nextlevel[inner_index].second);
					++inner_index;
				}
				new_inner->append((NodeBase *)(nextlevel[inner_index].second));
	
				nextlevel[i].first = nextlevel[inner_index].first;
				nextlevel[i].second = (uint64_t)new_inner;
				++inner_index;
				num_children -= bulk_num + 1;
			}
		}
		root = (NodeBase*)nextlevel[0].second;
		delete[] nextlevel;
	}

	void Bulkload(std::pair<Key, Value> * recs_, int num){
		uint64_t num_leaves = (num + BTreeLeaf<Key,Value>::maxEntries - 1) / BTreeLeaf<Key,Value>::maxEntries;
		uint64_t num_items = num;

		std::pair<Key, uint64_t> *nextlevel = new std::pair<Key, uint64_t>[num_leaves];
		
		
		BTreeLeaf<Key, Value>* last_leaf = nullptr;
		uint64_t cur = 0;
		for (size_t i = 0; i < num_leaves; ++i){
			BTreeLeaf<Key, Value>* leaf = new BTreeLeaf<Key, Value>();
			uint64_t bulk_num = num_items /(num_leaves - i);
			for(size_t j = 0; j < bulk_num; ++j){
				leaf->append(recs_[cur].first, recs_[cur].second);
				++cur;
			}
			nextlevel[i].first = recs_[cur-1].first; 
			nextlevel[i].second = (uint64_t)leaf;
			if(last_leaf != nullptr){
				last_leaf->next_leaf = leaf;	
			}
			else
				first_leaf = leaf;
			last_leaf = leaf;

			num_items -= bulk_num; 
		}
		
		if(first_leaf == last_leaf){
			root = first_leaf;
			return;
		}

		//build inner nodes
		// uint64_t num_children = num_leaves;
		uint64_t num_parents = num_leaves;
		while (num_parents != 1){
			uint64_t num_children = num_parents;
			num_parents = (num_children + BTreeInner<Key>::maxEntries - 1) / (BTreeInner<Key>::maxEntries - 1);
			uint64_t inner_index = 0;
			for (size_t i = 0; i < num_parents; i++){
				uint64_t bulk_num = num_children /(num_parents - i) - 1;
				BTreeInner<Key>* new_inner = new BTreeInner<Key>();
				for (size_t j = 0; j < bulk_num; j++){
					new_inner->append(nextlevel[inner_index].first, (NodeBase *)nextlevel[inner_index].second);
					++inner_index;
				}
				new_inner->append((NodeBase *)(nextlevel[inner_index].second));
				
				nextlevel[i].first = nextlevel[inner_index].first;
				nextlevel[i].second = (uint64_t)new_inner;
				++inner_index;
				num_children -= bulk_num + 1;
			}
		}
		root = (NodeBase*) nextlevel[0].second;
		delete[] nextlevel;

	}

    void makeRoot(Key k, NodeBase* leftChild,NodeBase* rightChild) {
        auto inner = new BTreeInner<Key>();
        inner->count = 1;
        inner->keys[0] = k;
        inner->children[0] = leftChild;
        inner->children[1] = rightChild;
        root = inner;
    }

    bool Insert(const Key &k, const Value &v) {

		// Current node
		NodeBase* node = root;

		// Parent of current node
		BTreeInner<Key>* parent = nullptr;

		while (node->type==PageType::BTreeInner) {
			auto inner = static_cast<BTreeInner<Key>*>(node);

			// Split eagerly if full
			if (inner->isFull()) {
			
			// Split
				Key sep; BTreeInner<Key>* newInner = inner->split(sep);
				
				if (parent)
					parent->insert(sep,newInner);
				else
					makeRoot(sep,inner,newInner);
				if(k >= sep)	inner = newInner;
			}
			parent = inner;
			node = inner->children[inner->lowerBound(k)];
			
		}

		auto leaf = static_cast<BTreeLeaf<Key,Value>*>(node);

      // Split leaf if full
      	if (leaf->count==leaf->maxEntries) {
        // Split
			Key sep; BTreeLeaf<Key,Value>* newLeaf = leaf->split(sep);
			if(k > sep) newLeaf->insert(k, v);
			else leaf->insert(k, v);
			if (parent)
				parent->insert(sep, newLeaf);
			else
				makeRoot(sep, leaf, newLeaf);
		} else {
			leaf->insert(k, v);
		}
		return true; // success
    }
	uint64_t FindInsert(const Key &key, Value val){
		return 0;
	}

    bool Find(const Key &k, Value& result) {
	
		#ifndef NDEBUG
        of_search_time++;
        #endif
		NodeBase* node = root;
		// Parent of current node
		BTreeInner<Key>* parent = nullptr;

		while (node->type==PageType::BTreeInner) {
			#ifndef NDEBUG
			of_probe_length++;
			#endif
			auto inner = static_cast<BTreeInner<Key>*>(node);
			node = inner->children[inner->lowerBound(k)];
		}

		BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
		unsigned pos = leaf->lowerBound(k);
		bool success;
		if ((pos < leaf->count) && (leaf->keys[pos]==k)) {
			success = true;
			result = leaf->payloads[pos];
		}
		return success;
    }

	bool Locate(const Key &k, Value* &result) {
	
		#ifndef NDEBUG
        of_search_time++;
        #endif
		NodeBase* node = root;
		// Parent of current node
		BTreeInner<Key>* parent = nullptr;

		while (node->type==PageType::BTreeInner) {
			#ifndef NDEBUG
			of_probe_length++;
			#endif
			auto inner = static_cast<BTreeInner<Key>*>(node);
			node = inner->children[inner->lowerBound(k)];
		}

		BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
		unsigned pos = leaf->lowerBound(k);
		bool success;
		if ((pos < leaf->count) && (leaf->keys[pos]==k)) {
			success = true;
			result = &leaf->payloads[pos];
		}
		return success;
    }

	bool Update(const Key &k, const Value& v) {

      NodeBase* node = root;
      // Parent of current node
      BTreeInner<Key>* parent = nullptr;

      while (node->type==PageType::BTreeInner) {
        auto inner = static_cast<BTreeInner<Key>*>(node);
        node = inner->children[inner->lowerBound(k)];
      }

      BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
      unsigned pos = leaf->lowerBound(k);
      bool success;
      if ((pos < leaf->count) && (leaf->keys[pos]==k)) {
        success = true;
        leaf->payloads[pos] = v;
      }
      return success;
    }

	bool Remove(const Key &k){
		return false;
	}

    uint64_t Scan(Key k, int range, Record<Key, Value>* output) {

      NodeBase *node = root;

    while (node->type == PageType::BTreeInner) {
		auto inner = static_cast<BTreeInner<Key> *>(node);
		node = inner->children[inner->lowerBound(k)];
    }

    auto leaf = static_cast<BTreeLeaf<Key, Value> *>(node);
    unsigned pos = leaf->lowerBound(k);
    int count = 0;

    while (leaf && count < range) {
      for (unsigned i = pos; i < leaf->count && count < range; i++) {
        	output[count++] = Record<Key, Value>(leaf->keys[i], leaf->payloads[i]);
      }

      if (count == range) {
        // scan() finishes at [leaf]
        break;
      } else {
        // proceed with next leaf
        auto next_leaf = leaf->next_leaf;
       
        if (!next_leaf) {
          // scan() finishes at [leaf]
          break;
        }
        leaf = next_leaf;
        pos = 0;
      }
    }
    return count;
    }

	uint64_t Size(){
		return 0;
	}

	inline Record<Key, Value> PopFront(){
        Record<Key, Value> result;
        
        return result;
    }
    inline Record<Key, Value> Front(){
        Record<Key, Value> result;
        
        return result;
    }
	bool IsEmpty(){
		return false;
	}

	bool IsFull(){
		return false;
	}

	void Prefetch(){
		return ;
	}

	void Dump(_key_t *keys, _value_t* values, uint64_t &start){
		return ;
	}

	uint64_t Memory(){
		return 0;
	}
  };


}//namespace btree
}//namespace lisk
#endif