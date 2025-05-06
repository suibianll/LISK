/*
 * @Author: Chu Zhaole
 * @Date: 2024-11-04 06:31:01
 * @LastEditTime: 2024-11-28 05:06:15
 * @LastEditors: Chu Zhaole
 * @Description: 
 * @FilePath: /LISK/src/alexnode.cc
 * Copyright (c) Chu Zhaole.
 */
#include "alexnode.h"
// #include "lisk.h"

namespace lisk{
ALEXNode::ALEXNode(_key_t *keys, _value_t *values, uint64_t num, uint16_t d){
    depth = d;
    prefix_len = 0;
    key_num = num;
    node_type = SUBALEX;
    capacity_ = num;
    idx = new(buffer) alex::Alex<_key_t, _value_t>();
    std::vector<std::pair<_key_t, _value_t>> loads;
    loads.reserve(num);
    for (size_t i = 0; i < num; i++){
        loads.emplace_back(keys[i], values[i]);
    }
    idx->bulk_load(loads.data(), num);
    
}
ALEXNode::ALEXNode(){
    idx = new(buffer) alex::Alex<_key_t, _value_t>();
}

bool ALEXNode::Insert(_key_t key, _value_t &val){
    auto ret = idx->insert(key, val);
    return ret.second;
}

bool ALEXNode::Find(_key_t key, _value_t &val){
    auto ret = idx->get_payload(key);
    if(ret){
        val = *ret;
        return true;
    }
    return false;
}
bool ALEXNode::Update(_key_t key, _value_t &val){
    auto ret = idx->update(key, val);
    return ret;
}
    
bool ALEXNode::Remove(_key_t key){
    return idx->erase_one(key);
}

bool ALEXNode::Scan(_key_t startKey, uint64_t len, _value_t* &result){
    return false;
}


}
