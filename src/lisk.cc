
#include "lisk.h"
namespace lisk
{
uint64_t index_depth = 0;
uint64_t inner_probe_length = 0;
uint64_t inner_insert = 0;
uint64_t leaf_probe_length = 0;
uint64_t of_probe_length = 0;
uint64_t of_search_time = 0;
uint64_t retrain_time = 0;
uint64_t hot_search = 0;
uint64_t btree_search = 0;
uint64_t learned_search = 0;
uint64_t verify_num = 0;
uint64_t total_leaf_num = 0;
uint64_t fmcd_D = 0;
uint64_t fmcd_success = 0;
uint64_t fmcd_failure = 0;
uint64_t expand_failure = 0;
uint64_t search_depth[1000]={0};
double latency_breakdown[100]={0};

uint64_t bytes[8][256];

Node* BTreeBulkloadTwo(Key &k1, _value_t v1, Key &k2, _value_t v2, int depth){

    _key_t slice_k1 = k1.getSlice(depth);
    _key_t slice_k2 = k2.getSlice(depth);
    BTreeLeafNode* new_layer = new BTreeLeafNode(depth);
    if(slice_k1 == slice_k2){//do not use path compression, for simplicity of subsquent insetions.
        Node* next_layer = BTreeBulkloadTwo(k1, v1, k2, v2, depth + SLICE);
        new_layer->append(slice_k1, (_value_t)next_layer);
    }
    else{
        
        if(slice_k1 < slice_k2){
            new_layer->append(slice_k1, v1);
            new_layer->append(slice_k2, v2);
        }
        else{
            new_layer->append(slice_k2, v2);
            new_layer->append(slice_k1, v1);
        }
    }
    return (Node *)(SET_TYPE(new_layer, BLEAF));
}

Node* BTreeBulkload(Key* keys, _value_t* values, int num, int depth){

    std::vector<_key_t> slice_ks;
    std::vector<_value_t> slice_vs;
        
    slice_ks.reserve(num);
    slice_vs.reserve(num);
    for (size_t i = 0; i < num; i++){
        _key_t k = keys[i].getSlice(depth);
        if(i < num - 1 && k == keys[i + 1].getSlice(depth)){
            size_t j = i + 1;
            while(j < num && keys[j].getSlice(depth) == k) ++j;
            Node* sub_index = Bulkload(&keys[i], &values[i], 0, j - i, depth + SLICE);//bulkload the subset
            slice_vs.push_back((_value_t)sub_index);//to fix, we should use the pointer to indicate the pointer is node or actural value
            i = j - 1;
        }
        else{
            slice_vs.push_back(values[i]);
        }
        slice_ks.push_back(k);
    }

    uint64_t slice_num = slice_ks.size();

    return BTree(slice_ks.data(), slice_vs.data(), slice_num, depth);
}

Node* ModelBulkLoad(Key* keys, _value_t* values, int num, int depth){
    int cpl = commonPrefixLen(keys[0], keys[num - 1]);

    std::vector<_key_t> slice_ks;
    std::vector<_value_t> slice_vs;
    
    slice_ks.reserve(num);
    slice_vs.reserve(num);
    for (size_t i = 0; i < num; i++){
        _key_t k = keys[i].getSlice(depth);

        if(i < num - 1 && k == keys[i + 1].getSlice(depth)){
            size_t j = i + 1;
            while(j < num && keys[j].getSlice(depth) == k) ++j;
            Node* sub_index = Bulkload(&keys[i], &values[i], 0, j - i, depth + SLICE);//bulkload the subset
            slice_vs.push_back((_value_t)sub_index);
            i = j - 1;
        }
        else{
            slice_vs.push_back(values[i]);
        }
        slice_ks.push_back(k);
    }

    uint64_t slice_num = slice_ks.size();
    if(slice_num < 2048){//too few unique keys, do not allow a single leaf node, this will complicate the process of insert operation
        return BTreeBulkload(keys, values, num, depth);
    }
    int seg_num = slice_num / 2000 > 0 ? slice_num / 2000 : 1;
    SecondDerivative<_key_t, _value_t> seg(slice_ks.data(), slice_num, seg_num);
    uint64_t leaf_num = seg.number();
    
    std::vector<_key_t> leaf_index;
    leaf_index.reserve(leaf_num);
    std::vector<_value_t> leaf_nodes; 
    leaf_nodes.reserve(leaf_num);

    uint64_t start_pos = 0;
    for (size_t i = 0; i < leaf_num; i++){
        int key_num = seg[i].number();
        leaf_nodes.emplace_back((_value_t)(new LeafNode(&slice_ks[start_pos], &slice_vs[start_pos], key_num, depth)));
        leaf_index.push_back(i == 0 ? std::numeric_limits<_key_t>::min() :slice_ks[start_pos]);
        start_pos += key_num;
    }
    for (size_t i = 0; i < leaf_nodes.size() - 1; i++){
        ((LeafNode *)leaf_nodes[i])->SetSibling(leaf_nodes[i + 1]);
    }
    
    Node* res = new InnerNode(leaf_index.data(), leaf_nodes.data(), leaf_nodes.size(), depth);

    return (Node *) (SET_TYPE(res, MINNER));
}

Node* Bulkload(Key* keys, _value_t* values, int left, int right, int depth){
    int num = right - left;
    if(num == 1){
        std::cout<<"bulkload error"<<std::endl;
        return nullptr;
    }
    else if(num < 4096){
        Node* node = BTreeBulkload(&keys[left], &values[left], num, depth);
        return node;
    }
    else{
        Node* node = ModelBulkLoad(&keys[left], &values[left], num, depth);
        return node;
    }
}

} // namespace lisk
