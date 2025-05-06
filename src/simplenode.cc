
#include "node.h"

namespace lisk{
SimpleLeafNode::SimpleLeafNode(Key *keys, _value_t *values, uint64_t num){
    assert(num >= 1 && num <= 16);
    int cpl = commonPrefixLen(keys[0], keys[num - 1]);
    depth = cpl / SLICE * SLICE;
    prefix_len = cpl - depth;
    key_num = num;
    node_type = SLEAF;
    capacity_ = 16;
    dataslot = (Key**) malloc(capacity_ * sizeof(Key *));
    for (size_t i = 0; i < num; i++){
        dataslot[i] = setHashKey(&keys[i]);
    }
    
}

SimpleLeafNode::SimpleLeafNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2){
    int cpl = commonPrefixLen(*(const_cast<Key*>(&k1)), *(const_cast<Key*>(&k2)));
    depth = cpl / SLICE * SLICE;
    prefix_len = cpl - depth;
    key_num = 2;
    node_type = SLEAF;
    capacity_ = 16;
    dataslot = (Key**) malloc(capacity_ * sizeof(Key *));
    if(k1 < k2){
        dataslot[0] = setHashKey(&k1);
        dataslot[1] = setHashKey(&k2);
    }
    else{
        dataslot[0] = setHashKey(&k2);
        dataslot[1] = setHashKey(&k1);
    }
}

bool SimpleLeafNode::Insert(Key &key, _value_t &val){
    int i = 0;
    for (; i < key_num; i++){
        Key* k = (Key *)GET_POINTER(dataslot[i]);
        if(*k > key) {
            break;
        }
    }
    if(i < key_num){
        if(unlikely(i == 0)){
            int cpl = commonPrefixLen(key, *(Key *)GET_POINTER(dataslot[key_num - 1]));
            depth = cpl / SLICE * SLICE;
            prefix_len = cpl - depth;
        }
        memmove(&dataslot[i + 1], &dataslot[i], sizeof(Key*) * (key_num -i));
        
    }
    else{
        int cpl = commonPrefixLen(key, *(Key *)GET_POINTER(dataslot[0]));
        depth = cpl / SLICE * SLICE;
        prefix_len = cpl - depth;
    }
    
    dataslot[i] = setHashKey(&key);
    key_num += 1;
        
    if(key_num == SLEAF_SIZE){//travert has bug
        std::vector<Key *> keys;
        std::vector<_value_t> values;
        keys.reserve(key_num);
        values.reserve(key_num);
        Dump(keys, values);
        ArtNode* new_inner = new ArtNode(keys.data(), values.data(), keys.size());
        SwapNodes(this, new_inner);

        delete ((SimpleLeafNode*)new_inner);
    }
    return true;
}

bool SimpleLeafNode::Find(Key &key, _value_t &val){
    uint16_t hash = hashKey(key);
    for (size_t i = 0; i < key_num; i++){
        if(hash != getHash(dataslot[i])) continue;
        Key* k = (Key *)GET_POINTER(dataslot[i]);
        if(key.verify(*k, depth)){
            val = k->val;
            return true;
        }
    }
    return false;
}
bool SimpleLeafNode::Update(Key &key, _value_t &val){
    uint16_t hash = hashKey(key);
    for (size_t i = 0; i < key_num; i++){
        if(hash != getHash(dataslot[i])) continue;
        Key* k = (Key *)GET_POINTER(dataslot[i]);
        if(*k == key){
            k->val=val;
            return true;
        }
    }
    return false;
}

bool SimpleLeafNode::Remove(Key &key){
    return false;
}

void SimpleLeafNode::Dump(std::vector<Key *> &keys, std::vector<_value_t> &values){
    for (size_t i = 0; i < key_num; i++){
        Key * k =(Key*)(GET_POINTER(dataslot[i]));
        keys.push_back(k);
        values.push_back(k->val);
    }
    
}
}
