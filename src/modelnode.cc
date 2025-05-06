/*
 * @Author: Chu Zhaole
 * @Date: 2024-11-25 07:06:31
 * @LastEditTime: 2025-05-06 08:31:21
 * @LastEditors: Chu Zhaole
 * @Description: 
 * @FilePath: /LISK/src/modelnode.cc
 * Copyright (c) Chu Zhaole.
 */
#include "modelnode.h"

namespace lisk{

bool ModelFind(Node* root, Key &key, _value_t &val, int depth){
    Node* n = root;
    _key_t k = key.getSlice(depth);
    uint64_t v;
    int search_depth = 0;
    while (n->node_type == MINNER){
        static_cast<InnerNode*>(n)->Find(k, v);
        n = reinterpret_cast<Node*>(GET_POINTER(v));
        #ifndef NDEBUG
        search_depth++;
        #endif
    }
    #ifndef NDEBUG
    if(depth == 0){
        index_depth += search_depth;
        learned_search ++;
    }
    #endif

    bool found = static_cast<LeafNode*>(n)->Find(k, val);

    return found;
}
Node* ModelInsert(Node* root, Key &key, _value_t val, int depth){
    Node* n = root;
    _key_t k = key.getSlice(depth);
    uint64_t v;

    while (n->node_type == MINNER){//forward to the leaf node
        static_cast<InnerNode*>(n)->Find(k, v);
        n = reinterpret_cast<Node*>(GET_POINTER(v));
    }
    std::pair<_key_t, _value_t> *split_leaf_index;
    uint64_t split_size = 0;
    bool splitIf = static_cast<LeafNode*>(n)->Insert(k, val, split_leaf_index, split_size);
    if(splitIf){
        for (int i = 0; i < split_size; i++){//if split
            #ifndef NDEBUG
                inner_insert ++;
            #endif
            reinterpret_cast<InnerNode*> (root)->Insert(split_leaf_index[i].first, split_leaf_index[i].second);
        }
    }
    return nullptr;
}

uint64_t ModelFindInsert(Node* root, Key &key, _value_t val, int depth){
    Node* n = root;
    Node* parent = n;
    _key_t k = key.getSlice(depth);
    uint64_t v;

    while (n->node_type == MINNER){//forward to the leaf node
        static_cast<InnerNode*>(n)->Find(k, v);
        parent = n;
        n = reinterpret_cast<Node*>(GET_POINTER(v));
    }
    std::pair<_key_t, _value_t> *split_leaf_index;
    uint64_t split_size = 0;
    uint64_t pos = static_cast<LeafNode*>(n)->FindInsert(k, val, split_leaf_index, split_size);
    if(!pos){//insert successfully
        for (int i = 0; i < split_size; i++){//if split
            reinterpret_cast<InnerNode*> (parent)->Insert(split_leaf_index[i].first, split_leaf_index[i].second);
        }
        return 0;
    }
    else{
        _value_t * p = reinterpret_cast<_value_t *>(pos);
        if(GET_TYPE(*p) == 0){
            Node* new_node = BTreeBulkloadTwo(*(reinterpret_cast<Key*>(GET_POINTER(*p))), (*p), (key), (val), depth + SLICE);
            *p = (uint64_t)new_node;
            return 0;
        }
    }
    return pos;
}

bool ModelUpdate(Node* root, Key &key, _value_t &val, int depth){
    Node* n = root;
    _key_t k = key.getSlice(depth);
    uint64_t v;

    while (n->node_type == MINNER){//forward to the leaf node
        static_cast<InnerNode*>(n)->Find(k, v);
        n = reinterpret_cast<Node*>(GET_POINTER(v));
    }
    bool success = static_cast<LeafNode*>(n)->Update(k, val);

    return success;
}
bool ModelRemove(Node* root, Key &key, int depth){
    return false;
}
uint64_t ModelScan(Key &key, int range, Record<_key_t, _value_t>* output, int depth){
    return 0;
}

}