#ifndef __LISK_MODELNODE_H__
#define __LISK_MODELNODE_H__
#include <cassert>
#include <cstring>
#include <atomic>
#include <immintrin.h>
#include <sched.h>
#include "node.h"
#include "util.h"

namespace lisk {

class InnerNode : public Node{
public:
    InnerNode(_key_t *keys, _value_t *values,  uint64_t num, uint16_t d = 0);
    InnerNode(Key *keys, _value_t *values,  uint64_t num);

    InnerNode();

    ~InnerNode(){
        delete[] dataslot;
        // free(dataslot);
    }

    bool Insert(_key_t &key, _value_t &val);
    
    bool Find(_key_t &key, _value_t &val);

    bool FindFirst(_value_t &val);

    bool Locate(_key_t &key, _value_t *&val);//locate the position of key, return the address of the value, used for the insert operation

    bool Update(_key_t &key, _value_t &val);

    bool Remove(_key_t &key);
    
    void Dump(std::vector<_key_t> &keys, std::vector<_value_t> &values);

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);
    uint64_t Size() const {return key_num;}

    uint64_t Memory() const;

    _key_t First(){return first_key_;}

constexpr static int BNODE_SIZE    = 12;
private:
    __always_inline int Predict(const _key_t & key) {
        double predict = (key - first_key_) * slope + intercept + 0.5;
        if(likely(predict >= 0 && predict < capacity_)) return predict;
        else if(predict >= capacity_)    return capacity_ - 1;
        else return 0;
        __builtin_unreachable(); 
    }
    double slope;
    double intercept;
    uint64_t capacity_;
    _key_t first_key_ = MAX_KEY;
    uint32_t of_num;
    uint8_t prefix[12];

private:
    Record<_key_t, _value_t> *dataslot;

}__attribute((packed));

class LeafNode : public Node{

public:
    LeafNode(Key *keys, _value_t *values, uint64_t num, uint16_t d = 0);
    LeafNode(_key_t *keys, _value_t *values, uint64_t num, uint16_t d = 0);
    LeafNode(_key_t *keys, _value_t *values, NodeParamer &p, uint16_t d = 0);

    LeafNode();
    ~LeafNode(){
        delete[] dataslot;
        // free(dataslot);
    }

    bool Insert(_key_t &key, _value_t &val, std::pair<_key_t, _value_t> *&split_leaf_index, uint64_t &split_size);

    uint64_t FindInsert(_key_t &key, _value_t &val, std::pair<_key_t, _value_t> *&split_leaf_index, uint64_t &split_size);
    //try to insert, if the key already exists, return the position of the key

    bool Find(_key_t &key, _value_t &val);

    bool Locate(_key_t &key, _value_t* &val);//locate the position of key, return the address of the value, used for the insert operation

    bool Update(_key_t &key, _value_t &val);
    
    bool Remove(_key_t &key);

    bool Scan(_key_t &startKey, uint64_t len, _value_t* &result);

    bool Append(_key_t &key, _value_t &val, uint64_t pos, bool &shouldSplit);

    bool DeAppend(_key_t &key, uint64_t pos, bool &shouldMerge);

    void DoSplit(std::pair<_key_t, _value_t> *&leaf_index, uint64_t &split_size);
    
    void Dump(_key_t *keys, _value_t* values, uint64_t &start);

    void SetSibling(uint64_t sibling){sibling_ = reinterpret_cast<LeafNode*>(sibling) ;}

    

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);

    uint64_t Size() const {return key_num;};

    uint64_t Memory() const;

private:
    
    __always_inline int Predict(const _key_t & key) {
        double predict = (key - first_key_) * slope + intercept + 0.5;
        if(likely(predict >= 0 && predict < capacity_)) return predict;
        else if(predict >= capacity_)    return capacity_ - 1;
        else return 0;
        __builtin_unreachable();
    }


public:
    double slope;
    double intercept;
    uint32_t capacity_;
    _key_t first_key_ = MAX_KEY;
    uint32_t of_num;
    uint8_t prefix[8];
    LeafNode * sibling_;//64 bytes
private:    

    Record<_key_t, _value_t> *dataslot;
 
}__attribute((packed));


inline void SwapNodes(InnerNode *oldone, InnerNode *newone){
    static const uint64_t HEADER_SIZE = sizeof(InnerNode);
    char * tmp = new char[HEADER_SIZE];
    sizeof(ArtNode);
    memcpy(tmp, oldone, HEADER_SIZE); // record the old node
    memcpy(oldone, newone, HEADER_SIZE); // replace the old node with the newone
    memcpy(newone, tmp, HEADER_SIZE);
}

inline void SwapNodes(LeafNode *oldone, LeafNode *newone){
    static const uint64_t HEADER_SIZE = sizeof(LeafNode);
    char * tmp = new char[HEADER_SIZE];

    memcpy(tmp, oldone, HEADER_SIZE); // record the old node
    memcpy(oldone, newone, HEADER_SIZE); // replace the old node with the newone
    memcpy(newone, tmp, HEADER_SIZE);
}

extern bool ModelFind(Node* root, Key &key, _value_t &val, int depth);
extern Node* ModelInsert(Node* root, Key &key, _value_t val, int depth);
extern uint64_t ModelFindInsert(Node* root, Key &key, _value_t val, int depth);
extern bool ModelUpdate(Node* root, Key &key, _value_t &val, int depth);
extern bool ModelRemove(Node* root, Key &key, int depth);
extern uint64_t ModelScan(Key &key, int range, Record<_key_t, _value_t>* output, int depth);

}

#endif