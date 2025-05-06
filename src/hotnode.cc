

#include "node.h"

namespace lisk{
HotNode::HotNode(Key *keys, _value_t *values, uint64_t num, uint16_t d){
    int cpl = commonPrefixLen(keys[0], keys[num - 1]);
    depth = d;
    prefix_len = cpl - depth;
    key_num = num;
    node_type = SUBHOT;
    capacity_ = num;
    idx = new(buffer) HOTIndex();
    for (size_t i = 0; i < num; i++){
        keys[i].setVal(values[i]);
        idx->insert(&keys[i]);
    }
    
}
HotNode::HotNode(Key **keys, _value_t *values, uint64_t num, uint16_t d){
    int cpl = commonPrefixLen(*(keys[0]), *(keys[num - 1]));
    depth = d;
    prefix_len = cpl - depth;
    key_num = num;
    node_type = SUBHOT;
    capacity_ = num;
    idx = new(buffer) HOTIndex();
    for (size_t i = 0; i < num; i++){
        keys[i]->setVal(values[i]);
        idx->insert(keys[i]);
    }
    
}

HotNode::HotNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2, uint16_t d){
    int cpl = commonPrefixLen(k1, k2);
    depth = d;
    prefix_len = cpl - depth;
    key_num = 2;
    node_type = SUBHOT;
    capacity_ = 2;
    idx = new(buffer) HOTIndex();
    idx->insert(&k1);
    idx->insert(&k2);
}
HotNode::HotNode(){
    idx = new(buffer) HOTIndex();
}

bool HotNode::Insert(Key &key, _value_t &val){
    return idx->insert(&key);
}

bool HotNode::Find(Key &key, _value_t &val){
    auto ret = idx->lookup(key.data);
    if(ret.mIsValid){
        val = ret.mValue->val;
        return true;
    }
    return false;
}
bool HotNode::Update(Key &key, _value_t &val){
    key.setVal(val);
    auto ret = idx->upsert(&key);
    if(ret.mIsValid) return true;
    return false;
}
    
bool HotNode::Remove(Key &key){
    return idx->remove(key.data);
}

bool HotNode::Scan(Key &startKey, uint64_t len, _value_t* &result){
    return false;
}
}
