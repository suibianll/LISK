/*
 * @Author: Chu Zhaole
 * @Date: 2024-01-12 08:34:45
 * @LastEditTime: 2025-05-06 08:37:42
 * @LastEditors: Chu Zhaole
 * @Description: 
 * @FilePath: /LISK/src/leafnode.cc
 * Copyright (c) Chu Zhaole.
 */

#include "modelnode.h"
// #include "lisk.h"

namespace lisk
{

uint64_t bulkload_error;
LeafNode::LeafNode(Key *keys, _value_t *values, uint64_t num, uint16_t d){
    assert(num > 1);
    int cpl = commonPrefixLen(keys[0], keys[num - 1]);
    std::vector<_key_t> slice_ks;
    std::vector<_value_t> slice_vs;
    slice_ks.reserve(num);
    slice_vs.reserve(num);
    for (size_t i = 0; i < num; i++){
        _key_t k = keys[i].getSlice(cpl);
        if(i < num - 1 && k == keys[i + 1].getSlice(cpl)){
            size_t j = i + 1;
            while(j < num && keys[j].getSlice(cpl) == k) ++j;
            Node* sub_index = Bulkload(&keys[i], &values[i], 0, j - i);//bulkload the subset
            slice_vs.push_back(SET_TYPE(sub_index, NEXTLEVEL));//to fix, we should use the pointer to indicate the pointer is node or actural value
            i = j - 1;
        }
        else{
            slice_vs.push_back(values[i]);
        }
        slice_ks.push_back(k);
    }
    depth = cpl;
    prefix_len = cpl;
    key_num = slice_ks.size();
    node_type = MLEAF;
    capacity_ = key_num * SCALE_FACTOR;
    first_key_ = slice_ks[0];
    
    dataslot = new Record<_key_t, _value_t>[capacity_];
    
    for (size_t i = 0; i < capacity_; i++){
        dataslot[i].key = FREE_FLAG;
        dataslot[i].val = 0;
    }
    FMCD<_key_t> model;
    model.train(slice_ks.data(), key_num);
    slope = model.a_ * SCALE_FACTOR;
    intercept = model.b_ * SCALE_FACTOR;

    int i = 0, last_i = i;
    uint64_t unique_key_num = slice_ks.size();
    int cur_pos = Predict(slice_ks[i]);
    for (i = 1; i < unique_key_num; i++){
        int predict = Predict(slice_ks[i]);
        if(predict != cur_pos){
            int count = i - last_i;
            if(count == 1){//no conflicts
                dataslot[cur_pos] = Record<_key_t, _value_t>(slice_ks[last_i], slice_vs[last_i]);
            }
            else{//conflicts
                dataslot[cur_pos].key = slice_ks[last_i];
                #ifndef NDEBUG
                if(count > SIMPLE_LEAF_NODE_LEVEL3_SIZE) 
                    expand_failure++;
                #endif
                OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&slice_ks[last_i], &slice_vs[last_i], count);
                dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);

            }
            last_i = i;
            cur_pos = predict;
        }
    }

    if(last_i != unique_key_num){//the rest part
        int count = unique_key_num - last_i;
        if(count == 1){//no conflicts
            dataslot[cur_pos] = Record<_key_t, _value_t>(slice_ks[last_i], slice_vs[last_i]);
        }
        else{//conflicts
            dataslot[cur_pos].key = slice_ks[last_i];
            OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&slice_ks[last_i], &slice_vs[last_i], count);
            dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
            
        }
    }
    
}

LeafNode::LeafNode(_key_t *keys, _value_t *values, NodeParamer &p, uint16_t d){

    depth = d;
    prefix_len = d;
    key_num = p.size;
    node_type = MLEAF;
    capacity_ = p.capacity;
    first_key_ = keys[0];
    of_num = 0;
     
    dataslot = new Record<_key_t, _value_t>[capacity_];

    for (size_t i = 0; i < capacity_; i++){
        dataslot[i].key = FREE_FLAG;
        dataslot[i].val = 0;
    }

    slope = p.a_;
    intercept = p.b_;

    int i = 0, last_i = i;
    int cur_pos = Predict(keys[i]);
    for (i = 1; i < p.size; i++){
        int predict = Predict(keys[i]);
        if(predict != cur_pos){
            int count = i - last_i;
            if(count == 1){//no conflicts
                dataslot[cur_pos] = Record<_key_t, _value_t>(keys[last_i], values[last_i]);
            }
            else{//conflicts
                dataslot[cur_pos].key = keys[last_i];
                OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&keys[last_i], &values[last_i], count);
                dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
                of_num += count;
            }
            last_i = i;
            cur_pos = predict;
        }
    }

    if(last_i != p.size){//the rest part
        int count = p.size - last_i;
        if(count == 1){//no conflicts
            dataslot[cur_pos] = Record<_key_t, _value_t>(keys[last_i], values[last_i]);
        }
        else{//conflicts
            dataslot[cur_pos].key = keys[last_i];
            OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&keys[last_i], &values[last_i], count);
            dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
            of_num += count;
        }
    }
}

LeafNode::LeafNode(_key_t *keys, _value_t *values, uint64_t num, uint16_t d){
    assert(num > 1);
    depth = d;
    prefix_len = d;
    key_num = num;
    node_type = MLEAF;
    capacity_ = key_num * SCALE_FACTOR;
    first_key_ = keys[0];
    
    dataslot = new Record<_key_t, _value_t>[capacity_];
    
    for (size_t i = 0; i < capacity_; i++){
        dataslot[i].key = FREE_FLAG;
        dataslot[i].val = 0;
    }

    LeafModel<_key_t> model;
    model.train(keys, key_num);
    slope = model.a_ * SCALE_FACTOR;
    intercept = model.b_ * SCALE_FACTOR;

    int i = 0, last_i = i;

    int cur_pos = Predict(keys[i]);
    for (i = 1; i < num; i++){

        int predict = Predict(keys[i]);
        if(predict != cur_pos){
            int count = i - last_i;
            if(count == 1){//no conflicts
                dataslot[cur_pos] = Record<_key_t, _value_t>(keys[last_i], values[last_i]);
            }
            else{//conflicts
                dataslot[cur_pos].key = keys[last_i];
                OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&keys[last_i], &values[last_i], count);
                dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
                of_num += count;
            }
            last_i = i;
            cur_pos = predict;
        }
    }

    if(last_i != num){//the rest part
        int count = num - last_i;
        if(count == 1){//no conflicts
            dataslot[cur_pos] = Record<_key_t, _value_t>(keys[last_i], values[last_i]);
        }
        else{//conflicts
            dataslot[cur_pos].key = keys[last_i];
            OFIndex<_key_t, _value_t> *new_node = new OFIndex<_key_t, _value_t>(&keys[last_i], &values[last_i], count);
            dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
            of_num += count;
        }
    }
}



LeafNode::LeafNode(){
    node_type = MLEAF;
    sibling_ = nullptr;
    return;
}


bool LeafNode::Append(_key_t &key, _value_t &val, uint64_t pos, bool &shouldSplit){
    if((dataslot[pos].key == FREE_FLAG) || (dataslot[pos].key == key)){
        dataslot[pos].key = key;
        dataslot[pos].val = val;
        return true;
    }
    else{
        bool success = false;
        if(GET_TYPE(dataslot[pos].val) != OFNODE){//no overflow before
            uint64_t initial_cap = 8;
            OFIndex<_key_t, _value_t> *of_tree = new OFIndex<_key_t, _value_t>(initial_cap);
            success = of_tree->Insert(dataslot[pos].key,dataslot[pos].val);
            success = of_tree->Insert(key,val);

            dataslot[pos].val = SET_TYPE(of_tree, OFNODE);
            of_num += 2;
            if(dataslot[pos].key > key)
                dataslot[pos].key = key;
        }
        else{//already overflow
            OFIndex<_key_t, _value_t> *of_tree = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[pos].val); 
            success = of_tree->Insert(key, val);
            if(dataslot[pos].key > key) dataslot[pos].key = key;
            of_num++;
            if((double)of_num / key_num > 0.9){
                shouldSplit=true;
            }
        }
        return success;
    }
    return true;
}


bool LeafNode::DeAppend(_key_t &key, uint64_t pos, bool &shouldMerge){
    if(dataslot[pos].key > key) return false;
    else{
        if(GET_TYPE(dataslot[pos].val) == VALUE){
            if(dataslot[pos].key == key){
                dataslot[pos].key = FREE_FLAG;
                return true;
            }
            return false;
        }
        else{
            OFIndex<_key_t, _value_t>* of_index = (OFIndex<_key_t, _value_t> *) GET_POINTER(dataslot[pos].val);
            bool succ = of_index->Remove(key);
            if(of_index->IsEmpty())
                dataslot[pos] = Record<_key_t, _value_t>();
            else if(dataslot[pos].key == key){
                dataslot[pos].key = of_index->Front().key;
            }
            return succ;
        }
    }
    return false;
}


bool LeafNode::Find(_key_t &key, _value_t &val){
    #ifndef NSTAT
    auto start = std::chrono::steady_clock::now();
    #endif
    
    if(unlikely(key < first_key_)){
        int predict =Predict(first_key_);
        if(unlikely(GET_TYPE(dataslot[predict].val) != OFNODE)) return false;
        else{
            OFIndex<_key_t, _value_t> * of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[predict].val);
            bool found =  of_index->Find(key, val);
            return found;
        }
    }

    int predict = Predict(key);

    if(unlikely(dataslot[predict].key > key)){
        return false;
    }
    else{
        if(unlikely(GET_TYPE(dataslot[predict].val) != OFNODE)){
            if((dataslot[predict].key == key)){
                val = dataslot[predict].val;
                #ifndef NSTAT
                auto end = std::chrono::steady_clock::now();
                latency_breakdown[LEAFPROBE] += (std::chrono::duration<double, std::micro>(end - start).count());
                #endif
                #ifndef NDEBUG
                leaf_probe_length++;
                #endif
                return true;
            }
            return false;
        }
        else{//search the overflow
            #ifndef NSTAT
            auto start = std::chrono::steady_clock::now();
            #endif
            OFIndex<_key_t, _value_t> * of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[predict].val);
            of_index->Prefetch();

            bool found =  of_index->Find(key, val);

            #ifndef NSTAT
                auto end = std::chrono::steady_clock::now();
                latency_breakdown[OVERFLOWPROBE] += (std::chrono::duration<double, std::micro>(end - start).count());
            #endif
            return found;
        }
    }
    return false;

}

bool LeafNode::Locate(_key_t &key, _value_t* &val){
        #ifndef NSTAT
    auto start = std::chrono::steady_clock::now();
    #endif
    if(unlikely(key < first_key_)){
        int predict =Predict(first_key_);
        if(unlikely(GET_TYPE(dataslot[predict].val) != OFNODE)) return false;
        else{
            OFIndex<_key_t, _value_t> * of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[predict].val);
            bool found =  of_index->Locate(key, val);
            return found;
        }
    }

    int predict = Predict(key);
    if(unlikely(dataslot[predict].key > key)){
        return false;
    }
    else{
        if(unlikely(GET_TYPE(dataslot[predict].val) != OFNODE)){
            if(dataslot[predict].key == key){
                val = &dataslot[predict].val;
                #ifndef NSTAT
                auto end = std::chrono::steady_clock::now();
                latency_breakdown[LEAFPROBE] += (std::chrono::duration<double, std::micro>(end - start).count());
                #endif
                #ifndef NDEBUG
                leaf_probe_length++;
                #endif
                return true;
            }
            return false;
        }
        else{//search the overflow
            #ifndef NSTAT
            auto start = std::chrono::steady_clock::now();
            #endif
            OFIndex<_key_t, _value_t> * of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[predict].val);
            bool found =  of_index->Locate(key, val);
            of_index->Prefetch();

            #ifndef NSTAT
                auto end = std::chrono::steady_clock::now();
                latency_breakdown[OVERFLOWPROBE] += (std::chrono::duration<double, std::micro>(end - start).count());
            #endif
            return found;
        }
    }
    return false;
}

bool LeafNode::Insert(_key_t &key, _value_t &val, std::pair<_key_t, _value_t> *&split_leaf_index, uint64_t &split_size){
    int predict;
    if(unlikely(key < first_key_)){
        predict = Predict(first_key_);
    }
    else 
        predict = Predict(key);
    bool shouldSplit = false;
    bool success = Append(key, val, predict, shouldSplit);
    key_num ++;
    if(shouldSplit){//should split the node
    #ifndef NDEBUG
        retrain_time ++;
    #endif
        DoSplit(split_leaf_index, split_size);
        return true;
    }
    return false;
}



uint64_t LeafNode::FindInsert(_key_t &key, _value_t &val, std::pair<_key_t, _value_t> *&split_leaf_index, uint64_t &split_size){
    int pos;
    if(unlikely(key < first_key_)){
        pos = Predict(first_key_);
    }
    else 
        pos = Predict(key);
    bool shouldSplit = false;
    if(dataslot[pos].key == FREE_FLAG){
        dataslot[pos].key = key;
        dataslot[pos].val = val;
        key_num ++;
    }
    else{
        bool success = false;
        if(GET_TYPE(dataslot[pos].val) != OFNODE){//no overflow before
            if(dataslot[pos].key == key){//find the key, return the pos
                return uint64_t(&dataslot[pos].val);
            }
            OFIndex<_key_t, _value_t> *of_tree = new OFIndex<_key_t, _value_t>(8);
            of_tree->Insert(dataslot[pos].key,dataslot[pos].val);
            of_tree->Insert(key,val);
            dataslot[pos].val = SET_TYPE(of_tree, OFNODE);
            if(dataslot[pos].key > key)
                dataslot[pos].key = key;
            key_num ++;
        }
        else{//already overflow
            OFIndex<_key_t, _value_t> *of_tree = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[pos].val); 
            of_tree->Prefetch();
            uint64_t res = of_tree->FindInsert(key, val);
            if(res != 0) return res;////find the key, return the pos
            key_num ++;
            if(dataslot[pos].key > key) dataslot[pos].key = key;
            if(of_tree->IsFull()){
                DoSplit(split_leaf_index, split_size);
            }
        }
    }
    return 0;
}

bool LeafNode::Update(_key_t &key, _value_t &val){
    int predict;
    if(unlikely(key < first_key_)){
        predict = Predict(first_key_);
    }
    else
        predict = Predict(key);
    if(unlikely(dataslot[predict].key > key)){
        return false;
    }
    else{
        if(GET_TYPE(dataslot[predict].val) != OFNODE){
            if(dataslot[predict].key == key){
                dataslot[predict].val = val;
                return true;
            }
            return false;
        }
        else{//the overflow
            OFIndex<_key_t, _value_t> * of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[predict].val);
            of_index->Prefetch();
            return of_index->Update(key, val);
        }
    }
    return false;
}


bool LeafNode::Remove(_key_t &key){
    int predict = Predict(key);
    bool shoulMerge = false;
    DeAppend(key, predict, shoulMerge);
    key_num --;
    //todo: merge the node
    return true;
}


bool LeafNode::Scan(_key_t &startKey, uint64_t len, _value_t* &result){
    int predict = Predict(startKey);
    //todo: 
    return true;
}


void LeafNode::DoSplit(std::pair<_key_t, _value_t> *&split_leaf_index, uint64_t &split_size){
    _key_t* keys = new _key_t[key_num];
    _key_t* values = new _value_t[key_num];

    uint64_t start = 0;
    Dump(keys, values, start);


    uint64_t max_sd = 0;
     int no_seg_range = key_num / 4;
    uint64_t split_pos = no_seg_range;
    for (size_t i = no_seg_range; i < key_num - no_seg_range; i++){
        uint64_t sd = labs((keys[i] - keys[i-1]) - (keys[i-1] - keys[i-2]));
        if(sd > max_sd){
            max_sd = sd;
            split_pos = i - 1;
        }
    }
    int leaf_num = 2;
    LeafNode** leaf_nodes = new LeafNode*[leaf_num];
    split_leaf_index = new std::pair<_key_t,_value_t>[leaf_num - 1];
    leaf_nodes[0] = new LeafNode(&keys[0], &values[0], split_pos, depth);
    leaf_nodes[1] = new LeafNode(&keys[split_pos], &values[split_pos], key_num - split_pos, depth);
    split_leaf_index[0] = std::make_pair(keys[split_pos], (uint64_t)leaf_nodes[1]);

    split_size = leaf_num - 1;
    for (size_t i = 0; i < leaf_num - 1; i++){
        leaf_nodes[i]->sibling_ = leaf_nodes[i + 1];
    }
    leaf_nodes[leaf_num - 1]->sibling_ = sibling_;
    
    SwapNodes(this, leaf_nodes[0]);
    delete leaf_nodes[0];
}

void LeafNode::Dump(_key_t *keys, _value_t* values, uint64_t &start){
    for (size_t i = 0; i < capacity_; i++){
        if(GET_TYPE(dataslot[i].val)==OFNODE){
            ((OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[i].val))->Dump(keys, values, start);
        }
        else if(dataslot[i].key != FREE_FLAG){
            keys[start] = dataslot[i].key;
            values[start] = dataslot[i].val;
            start ++;
        }
    }
}
void LeafNode::Stat(uint64_t &key_num, uint64_t &overflow_key_num, double &fill_rate, double &overflow_rate, uint64_t &overflow_num){
    for (size_t i = 0; i < capacity_; i++ ){
        if( (GET_TYPE(dataslot[i].val)==OFNODE)){
            OFIndex<_key_t, _value_t> *of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[i].val);
            overflow_key_num += of_index->Size();
            ++overflow_num;
        }
        
    }
    key_num = this->key_num - overflow_key_num;

    fill_rate = (double)(this->key_num - overflow_key_num) / this->key_num;
    overflow_rate = (double)overflow_key_num / this->key_num;
}


uint64_t LeafNode::Memory() const{
    uint64_t res = sizeof(LeafNode);
    res += (capacity_ * sizeof(Record<_key_t, _value_t>));
    for (size_t i = 0; i < capacity_; i ++ ){
        if(GET_TYPE(dataslot[i].val) == OFNODE){
            OFIndex<_key_t, _value_t> *of_index = (OFIndex<_key_t, _value_t> *)GET_POINTER(dataslot[i].val);
            res += of_index->Memory();
        }
    }
    return res;
}

    
} // namespace lisk
