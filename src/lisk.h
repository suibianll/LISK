#ifndef __LISK_H__
#define __LISK_H__

#include <iostream>
#include <vector>
#include "modelnode.h"
#include "btreenode.h"
#include "node.h"

namespace lisk
{
#define GET_RADIX(x) (((uint64_t)x & 0xfff0000000000000ul) >> 52)
#define GET_KEY(x) (((uint64_t)x & (~0xfff0000000000000)))
#define GET(x,k) (((uint64_t)x & (0xff00000000000000ul>>(k*4*2)))>>(64ul-k*4*2-8))


class LISK
{
public:
    LISK();

    LISK(std::vector<Key> &keys, std::vector<uint64_t> &values, uint64_t number);

    ~LISK();

    bool find(Key &key, _value_t &val);
    
    bool insert(Key &key, _value_t &val);
    
    bool update(Key &key, _value_t &val);

    bool remove(Key &key);

    void bulk_load(std::vector<Key> &keys, std::vector<_value_t> &values, uint64_t number);

    void Stat();

private:

    Node *root_;
};


LISK::LISK(){
    root_ = nullptr;
}

LISK::LISK(std::vector<Key> &keys, std::vector<uint64_t> &values, uint64_t number)
{
    root_ = Bulkload(keys.data(), values.data(), 0, number);
}
void LISK::bulk_load(std::vector<Key> &keys, std::vector<_value_t> &values, uint64_t number){
    root_ = Bulkload(keys.data(), values.data(), 0, number);
}

LISK::~LISK()
{

}    


bool LISK::find(Key &key, _value_t &val){
    Node* n = (Node*)(GET_POINTER(root_));
    uint64_t v;
    int depth = 0;

    bool found = ModelFind(n, key, v, depth); 
    while (found && GET_TYPE(v) != 0){
        depth += SLICE;
        n = reinterpret_cast<Node*>(GET_POINTER(v));
        n->Prefetch();

        if(GET_TYPE(v) <= MLEAF)
            found = ModelFind(n, key, v, depth);
        else
            found = BTreeFind(n, key, v, depth);

    }
    Key* exist = reinterpret_cast<Key*>(GET_POINTER(v));

    if(found){
        if(depth + SLICE < key.len && !key.verify(*(exist), depth + SLICE)) return false;
        val = v;
        return true;
    }
    return false;
}

bool LISK::insert(Key &key, _value_t &val){
    Node* n = (Node*)(GET_POINTER(root_));
    Node* parent = n;
    uint64_t v;
    int depth = 0;

    bool found = ModelFind(n, key, v, depth); 
    while (found && GET_TYPE(v) != 0){
        depth += SLICE;
        parent = n;
        n = reinterpret_cast<Node*>(GET_POINTER(v));
        n->Prefetch();
        if(GET_TYPE(v) < BINNER)
            found = ModelFind(n, key, v, depth);
        else 
            found = BTreeFind(n, key, v, depth);
    }

    if(!found){//
        if(n->node_type < BINNER){//model insert
            ModelInsert(n, key, val, depth);
        }
        else{//Btree insert
            Node* new_root = BTreeInsert(n, key, val, depth);
            if(new_root){//insertion leads to root split
                _value_t new_val = SET_TYPE(new_root, BINNER);
                if(parent->node_type < BINNER) ModelUpdate(parent, key, new_val, depth - SLICE);
                else BTreeUpdate(parent, key, new_val, depth - SLICE);
            }
        }
        return true;
    }
    else{
        if(!key.verify(*(reinterpret_cast<Key*>(GET_POINTER(v))), depth + SLICE)){//found the slice key, but the true keys do not equal
            Node* new_node = BTreeBulkloadTwo(key, val, *(reinterpret_cast<Key*>(GET_POINTER(v))), v, depth + SLICE);//create sub-tree
            _value_t new_val = SET_TYPE(new_node, BLEAF);
            if(n->node_type < BINNER) ModelUpdate(n, key, new_val, depth);//update the val in this layer
            else BTreeUpdate(n, key, new_val, depth);
            return true;
        }
    }
    //the key has already existed
    return false;
}


bool LISK::update(Key &key, _value_t &val){
    Node* n = (Node*)(GET_POINTER(root_));
    uint64_t v;
    int depth = 0;
    bool found = ModelFind(n, key, v, depth); 
    while (found && GET_TYPE(v) != 0){
        depth += SLICE;
        n = reinterpret_cast<Node*>(GET_POINTER(v));
        n->Prefetch();
        if(GET_TYPE(v) < BINNER)
            found = ModelFind(n, key, v, depth);
        else 
            found = BTreeFind(n, key, v, depth);
    }
    Key* exist = reinterpret_cast<Key*>(GET_POINTER(v));
    if(found && key.verify(*(exist), depth + SLICE)){
        exist->setVal(val);
        return true;
    }
    return false;
}


bool LISK::remove(Key &key){
    return false;
}


void LISK::Stat(){
    uint64_t total_leaf_num = 0;
    double average_key_num=0, average_overflow_key_num=0, average_overflow_num=0;
    double average_fill_rate = 0, average_overflow_rate = 0;

    uint64_t key_num=0, overflow_key_num=0, overflow_num=0;
    double fill_rate, overflow_rate;

    uint64_t inner_mem = 0;
    uint64_t leaf_mem = 0;

    Node* n = (Node*)(GET_POINTER(root_));
    (reinterpret_cast<InnerNode*> (n))->Stat(key_num, overflow_key_num, fill_rate, overflow_rate, overflow_num);
    inner_mem = (reinterpret_cast<InnerNode*> (n))->Memory();
    std::cout<<"----------------------------root----------------------------"<<std::endl;
    std::cout<<"keys stored in dataslot:"<<key_num<<std::endl;
    std::cout<<"keys stored in overflow index:"<<overflow_key_num<<std::endl;
    std::cout<<"fill rate:"<<fill_rate<<std::endl;
    std::cout<<"overflow rate:"<<overflow_rate<<std::endl;
    std::cout<<"overflow index number:"<<overflow_num<<std::endl;
    std::cout<<"average key number in overflow index:"<<(double) overflow_key_num / overflow_num<<std::endl;

    uint64_t v;
    uint64_t first;
    while (n->node_type == MINNER){
        InnerNode* inner = static_cast<InnerNode*>(n);
        inner->FindFirst(v);
        n = reinterpret_cast<Node*>(GET_POINTER(v));
    }
    LeafNode* leaf = static_cast<LeafNode*>(n);
    average_key_num=0, average_overflow_key_num=0, average_overflow_num=0;
    average_fill_rate = 0, average_overflow_rate = 0;
    while(leaf){
        overflow_key_num=0, overflow_num=0;
        total_leaf_num ++;
        leaf_mem += leaf->Memory();
        leaf->Stat(key_num, overflow_key_num, fill_rate, overflow_rate, overflow_num);
        average_key_num += key_num;
        average_overflow_key_num += overflow_key_num;
        average_overflow_num += overflow_num;
        average_fill_rate += fill_rate;
        average_overflow_rate += overflow_rate;
        leaf=leaf->sibling_;
    }

    std::cout<<"----------------------------leaf----------------------------"<<std::endl;
    std::cout<<"leaf num:"<< total_leaf_num<<std::endl;
    std::cout<<"keys stored in dataslot:"<<average_key_num/total_leaf_num<<std::endl;
    std::cout<<"keys stored in overflow index:"<<average_overflow_key_num/total_leaf_num<<std::endl;
    std::cout<<"fill rate:"<<average_fill_rate/total_leaf_num<<std::endl;
    std::cout<<"overflow rate:"<<average_overflow_rate/total_leaf_num<<std::endl;
    std::cout<<"overflow index number:"<<average_overflow_num/total_leaf_num<<std::endl;
    std::cout<<"average key number in overflow index:"<<(double) average_overflow_key_num / average_overflow_num<<std::endl;

    std::cout<<"----------------------------Memory----------------------------"<<std::endl;
    std::cout<<"inner node size:"<< (double)inner_mem / 1024.0 /1024.0 << " MB"<<std::endl;
    std::cout<<"leaf node size:"<< (double)leaf_mem / 1024.0 /1024.0 << " MB"<<std::endl;
}

} // namespace lisk
#endif