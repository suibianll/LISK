
#include "node.h"

namespace lisk {
ArtNode::ArtNode(Key *keys, _value_t *values, uint64_t num, uint16_t d){
    idx.root = nullptr;
    idx.size = 0;
    art_tree_init(&idx);
    
    int cpl = commonPrefixLen(keys[0], keys[num - 1]);
    depth = d;
    prefix_len = cpl - depth;
    key_num = num;
    node_type = SUBART;
    capacity_ = num;
    for (size_t i = 0; i < num; i++){
        art_insert(&idx, keys[i].data + depth, keys[i].len - depth, (void *)values[i]);
    }
    
}
ArtNode::ArtNode(Key **keys, _value_t *values, uint64_t num, uint16_t d){
    idx.root = nullptr;
    idx.size = 0;

    int cpl = commonPrefixLen(*keys[0], *keys[num - 1]);
    depth = d;
    prefix_len = cpl - depth;
    key_num = num;
    node_type = SUBART;
    capacity_ = num;
    for (size_t i = 0; i < num; i++){
        art_insert(&idx, keys[i]->data + depth, keys[i]->len - depth, (void *)values[i]);
    }    
}

ArtNode::ArtNode(Key &k1, _value_t &v1, Key &k2, _value_t &v2, uint16_t d){
    idx.root = nullptr;
    idx.size = 0;
    int cpl = commonPrefixLen(k1, k2);
    depth = d;
    prefix_len = d;
    key_num = 2;
    node_type = SUBART;
    capacity_ = 2;
    art_insert(&idx, k1.data + depth, k1.len - depth, (void *)v1);
    art_insert(&idx, k2.data + depth, k2.len - depth, (void *)v2);
}
ArtNode::ArtNode(){
    idx.root = nullptr;
    idx.size = 0;
}

bool ArtNode::Insert(Key &key, _value_t &val){

    void * res = art_insert(&idx, key.data + depth, key.len - depth, (void *)val);
    return !res;
}

bool ArtNode::Find(Key &key, _value_t &val){
    val = (_value_t)art_search(&idx, key.data + depth, key.len - depth);
    return val != 0;
}

bool ArtNode::Update(Key &key, _value_t &val){
    void * res = art_insert(&idx, key.data + depth, key.len - depth, (void *)val);
    return !res;
}
    
bool ArtNode::Remove(Key &key){
    void * res = art_delete(&idx, key.data + depth, key.len - depth);
    return !res;
}

bool ArtNode::Scan(Key &startKey, uint64_t len, _value_t* &result){
    return false;
}


}