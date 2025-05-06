
#ifndef __LISK_NODE_H__
#define __LISK_NODE_H__
#include "util.h"
#include "ofnode.h"
#include "art.h"
#include "clht.h"
#include "btree.h"
// #include "modelnode.h"
// #include "btreenode.h"

#include "hot_src/HOTSingleThreaded.hpp"
#include "hot_src/HOTSingleThreadedInterface.hpp"
#include "alex_src/alex.h"

#include <cstring>
#include <iostream>
#include <vector>
namespace lisk{
#if defined(CLHT) 
template<class T, class P>
using OFIndex = clht<T,P>;
#elif defined(HASHMAP)
template<class T, class P>
using OFIndex = HashNode<T,P>;
#elif defined(ARRAY)
template<class T, class P>
using OFIndex = OFArray<T,P>;
#elif defined(BTREE)
template<class T, class P>
using OFIndex = btree::BTree<T,P>;
#elif defined(FLAT)
template<class T, class P>
using OFIndex = FlatOFNode<T,P>;
#endif



extern uint64_t bulkload_error;
class Node
{   
public:
    Node() = default;

    void SetDepth(int d){depth = d;}
    __always_inline void Prefetch(){prefetch(this);}
public:
    uint8_t depth;//node depth
    uint8_t node_type;
    uint16_t prefix_len;
    uint32_t key_num;
};



constexpr static int SLEAF_SIZE    = 16;
class SimpleLeafNode : public Node{

public:
    SimpleLeafNode(Key *keys, _value_t *values, uint64_t num);

    SimpleLeafNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2);
    ~SimpleLeafNode(){free(dataslot);}

    bool Insert(Key &key, _value_t &val);

    bool Find(Key &key, _value_t &val);

    bool Update(Key &key, _value_t &val);
    
    bool Remove(Key &key);

    bool Scan(Key &startKey, uint64_t len, _value_t* &result);

    bool Append(Key &key, _value_t &val, uint64_t pos, bool &shouldSplit);

    bool DeAppend(Key &key, uint64_t pos, bool &shouldMerge);
    
    void Dump(std::vector<Key *> &keys, std::vector<_value_t> &values);

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);

    uint64_t Size() const {return key_num;};

    uint64_t Memory() const {return key_num;}

public:
    uint32_t capacity_;
    uint8_t prefix[12];
private:    

    Key **dataslot;
 
}__attribute((packed));


class HotNode : public Node{

// Key extractor used in HOT
template <typename T> class KeyExtracter {
  public:
    using KeyType = uint8_t const *;
    __always_inline KeyType operator()(const T &value) { return (uint8_t const *) value->data; }
    __always_inline KeyType operator()(const KeyType k) { return k; }
};

using HOTIter = hot::singlethreaded::HOTSingleThreadedIterator<Key*>;
using HOTIndex = hot::singlethreaded::HOTSingleThreaded<Key*, KeyExtracter>;

public:
    HotNode(Key *keys, _value_t *values, uint64_t num, uint16_t d = 0);
    HotNode(Key **keys, _value_t *values, uint64_t num, uint16_t d = 0);
    HotNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2, uint16_t d = 0);

    HotNode();
    ~HotNode(){delete idx;}

    bool Insert(Key &key, _value_t &val);

    bool Find(Key &key, _value_t &val);

    bool Update(Key &key, _value_t &val);
    
    bool Remove(Key &key);

    bool Scan(Key &startKey, uint64_t len, _value_t* &result);

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);

    uint64_t Size() const {return 0;}

    uint64_t Memory() const {return idx->getStatistics().first;}

private:    
    uint32_t capacity_;
    uint8_t buffer[12];
    HOTIndex* idx;
 
}__attribute((packed));


class ArtNode : public Node{

public:
    ArtNode(Key *keys, _value_t *values, uint64_t num, uint16_t d = 0);
    ArtNode(Key **keys, _value_t *values, uint64_t num, uint16_t d = 0);
    ArtNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2, uint16_t d = 0);

    ArtNode();
    ~ArtNode(){art_tree_destroy(&idx);}

    bool Insert(Key &key, _value_t &val);

    bool Find(Key &key, _value_t &val);

    bool Update(Key &key, _value_t &val);
    
    bool Remove(Key &key);

    bool Scan(Key &startKey, uint64_t len, _value_t* &result);

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);

    uint64_t Size() const {return 0;}

    uint64_t Memory() const {return 0;}

private:    
    uint32_t capacity_;
    uint8_t prefix[4];
    art_tree idx;
 
};//__attribute((packed));



inline void SwapNodes(SimpleLeafNode *oldone, HotNode *newone){
    static const uint64_t HEADER_SIZE = sizeof(SimpleLeafNode);
    char * tmp = new char[HEADER_SIZE];

    memcpy(tmp, oldone, HEADER_SIZE); // record the old node
    memcpy(oldone, newone, HEADER_SIZE); // replace the old node with the newone
    memcpy(newone, tmp, HEADER_SIZE);
}

inline void SwapNodes(SimpleLeafNode *oldone, ArtNode *newone){
    static const uint64_t HEADER_SIZE = sizeof(SimpleLeafNode);
    char * tmp = new char[HEADER_SIZE];

    memcpy(tmp, oldone, HEADER_SIZE); // record the old node
    memcpy(oldone, newone, HEADER_SIZE); // replace the old node with the newone
    memcpy(newone, tmp, HEADER_SIZE);
}

extern Node* BTreeBulkloadTwo(Key &k1, _value_t v1, Key &k2, _value_t v2, int depth);
extern Node* BTreeBulkload(Key* keys, _value_t* values, int num, int depth);
extern Node* ModelBulkLoad(Key* keys, _value_t* values, int num, int depth = 0);
extern Node* Bulkload(Key* keys, _value_t* values, int left, int right, int depth = 0);

}

#endif