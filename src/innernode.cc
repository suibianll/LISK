
#include "modelnode.h"

namespace lisk
{
InnerNode::InnerNode(_key_t *keys, _value_t *values,  uint64_t num, uint16_t d){
    node_type = MINNER;
    first_key_ = keys[0];
    capacity_ = (num * SCALE_FACTOR * 4);
    key_num = num;
    depth = d;
    prefix_len = d;
    
    if(num < BNODE_SIZE){
        capacity_ = BNODE_SIZE;
        dataslot = new Record<_key_t, _value_t>[capacity_];

        key_num = num;
        for (size_t i = 0; i < num; i++){
            dataslot[i].key = keys[i];
            dataslot[i].val = values[i];
        }
        for (size_t i = num; i < capacity_; i++){
            dataslot[i].key = FREE_FLAG;
            dataslot[i].val = 0;
        }
        return;
    }
    dataslot = new Record<_key_t, _value_t>[capacity_];

    for (size_t i = 0; i < capacity_; i++){
        dataslot[i].key = FREE_FLAG;
        dataslot[i].val = 0;
    }

    InnerModel<_key_t> model;
    model.train(keys, num);

    slope = model.a_ * SCALE_FACTOR * 4;
    intercept = model.b_ * SCALE_FACTOR * 4;

    int i = 0, last_i = i;
    int cur_pos = Predict(keys[i]);
    for (i = 1; i < num; i++){
        int predict = Predict(keys[i]);
        if(predict != cur_pos){
            int count = i - last_i;

            if(count == 1){
                dataslot[cur_pos] = Record<_key_t,_value_t>(keys[last_i], values[last_i]);
            }
            else{
                dataslot[cur_pos].key = keys[last_i];
                InnerNode *new_node = new InnerNode(&keys[last_i], &values[last_i], count, d);
                dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
            }

            last_i = i;
            cur_pos = predict;
        }
    }
    if(last_i != num){//the rest part
        int count = num - last_i;
        if(count == 1){
            dataslot[cur_pos] = Record<_key_t,_value_t>(keys[last_i], values[last_i]);
        }
        else{
            dataslot[cur_pos].key = keys[last_i];
            InnerNode *new_node = new InnerNode(&keys[last_i], &values[last_i], count, depth);
            dataslot[cur_pos].val = SET_TYPE(new_node, OFNODE);
        }
    }
    //duplicate the slot to make the prediction precise, avoiding the local search in the inner node
    int start=0, end=0;
    while(dataslot[start].key==FREE_FLAG) start++;
    end = start + 1;
    while (end < capacity_){
        while(end < capacity_ && dataslot[end].key == FREE_FLAG) end ++;
        _key_t key = dataslot[start].key;
        _value_t val = SET_TYPE(dataslot[start].val, DUPVALUE);//it is ok because when the original type is OFNODE, it can also be converted to DUPOFNODE
        for (size_t i = start + 1; i < end; i++){
            dataslot[i].key = key;
            dataslot[i].val = val;
        }
        start = end;
        end++;
    }
}

bool InnerNode::Find(_key_t &key, _value_t &val){
    if(key_num < BNODE_SIZE){
        int i = 0;
        for (; i < key_num; i++){
            if(dataslot[i].key > key) break;
            #ifndef NDEBUG
            inner_probe_length++;
            #endif
        }
        val = dataslot[i - 1].val;
        return true;
    }
    int predict = Predict(key);
    while (unlikely(key < dataslot[predict].key)){
        predict --;
        #ifndef NDEBUG
        inner_probe_length++;
        #endif
    }
    val = dataslot[predict].val;
    return true;
}

bool InnerNode::FindFirst(_value_t &val){
    for (size_t i = 0; i < capacity_; i++){
        if(dataslot[i].key != FREE_FLAG){
            val = dataslot[i].val;
            return true;
        }
    }
    return false;
    
}

bool InnerNode::Locate(_key_t &key, _value_t* &val){
    if(key_num < BNODE_SIZE){
        int i = 0;
        for (; i < key_num; i++){
            if(dataslot[i].key > key) break;
            #ifndef NDEBUG
            inner_probe_length++;
            #endif
        }
        if(unlikely(i==0)) i=1;//out of range insert
        val = &(dataslot[i - 1].val);
        return true;
    }
    int predict;
    if(unlikely(key < first_key_)) predict = Predict(first_key_);//out of range insert
    else predict = Predict(key);
    while ((key < dataslot[predict].key)){
        predict --;
        #ifndef NDEBUG
        inner_probe_length++;
        #endif
    }
    val = &(dataslot[predict].val);
    return true;
}


bool InnerNode::Insert(_key_t &key, _value_t &val){
    if(key_num < BNODE_SIZE){
        int i = 0;
        for (; i < key_num; i++){
            if(dataslot[i].key > key) break;
        }
        memmove(&dataslot[i + 1], &dataslot[i], sizeof(Record<_key_t,_value_t>) * (key_num -i));
        dataslot[i] = Record<_key_t,_value_t>(key, val);
        key_num += 1;
        
        if(key_num == BNODE_SIZE){
            std::vector<_key_t> keys;
            std::vector<_value_t> values;
            keys.reserve(key_num);
            values.reserve(key_num);
            Dump(keys, values);
            InnerNode *new_inner = new InnerNode(keys.data(), values.data(), keys.size(), depth);
            SwapNodes(this, new_inner);
            delete new_inner;
        }
        return true;
    }
    int predict = Predict(key);
    
    if(dataslot[predict].key == FREE_FLAG){
        dataslot[predict] = Record<_key_t,_value_t>(key,val);
        int end = predict + 1;
        _key_t kk = key;
        _value_t vv = SET_TYPE(val, DUPVALUE);
        while (end < capacity_ && dataslot[end].key == FREE_FLAG){
            dataslot[end] = Record<_key_t,_value_t>(kk, vv);
            end ++;
        }
        key_num += 1;
        return true;
    }
    else{
        while ((key < dataslot[predict].key)){
            predict --;
            #ifndef NDEBUG
            inner_probe_length++;
            #endif
        }
        int type = GET_TYPE(dataslot[predict].val);
        if(type == VALUE){//point to leaf
            _key_t keys[2];
            _value_t values[2];
            if(dataslot[predict].key > key){
                keys[0] = key;
                values[0] = val;
                keys[1] = dataslot[predict].key;
                values[1] = dataslot[predict].val;
            }
            else{
                keys[1] = key;
                values[1] = val;
                keys[0] = dataslot[predict].key;
                values[0] = dataslot[predict].val;
            }

            InnerNode * new_inner = new InnerNode(keys, values, 2);
            dataslot[predict].key = keys[0];
            dataslot[predict].val = SET_TYPE(new_inner, OFNODE);
            int end = predict + 1;
            while (end < capacity_ && GET_TYPE(dataslot[end].val) == DUPVALUE){
                dataslot[end] = Record<_key_t,_value_t>(keys[0], SET_TYPE(new_inner, DUPOFNODE));
                end ++;
            }
            key_num += 1;
            return true;
        }
        else if(type == OFNODE){
            InnerNode * inner = (InnerNode *)GET_POINTER(dataslot[predict].val);
            key_num += 1;
            return inner->Insert(key,val);
        }
        else if(type == DUPVALUE|| type == DUPOFNODE){
            dataslot[predict] = Record<_key_t,_value_t>(key,val);
            int end = predict + 1;
            _key_t kk = key;
            _value_t vv = SET_TYPE(val, DUPVALUE);
            while (end < capacity_ && (GET_TYPE(dataslot[end].val) == DUPVALUE||GET_TYPE(dataslot[end].val) == DUPOFNODE)){
                dataslot[end] = Record<_key_t,_value_t>(kk, vv);
                end ++;
            }
            key_num += 1;
            return true;
        }
    }
    return false;
    
}


void InnerNode::Dump(std::vector<_key_t> &keys, std::vector<_value_t> &values){
    for (size_t i = 0; i < capacity_; i++ ){
        if((GET_TYPE(dataslot[i].val)==OFNODE)){
            Node* n = (Node*) GET_POINTER(dataslot[i].val);
            if(n->node_type == MINNER){
                InnerNode* inner = (InnerNode*)n;
                inner->Dump(keys, values);
            }
            
        }
        else if(GET_TYPE(dataslot[i].val)==VALUE){
            keys.push_back(dataslot[i].key);
            values.push_back(dataslot[i].val);
        }
    }
}


void InnerNode::Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num){
    for (size_t i = 0; i < capacity_; i++ ){
        if((GET_TYPE(dataslot[i].val)==OFNODE)){
            Node* n = (Node*) GET_POINTER(dataslot[i].val);
            if(n->node_type == MINNER){
                ++overflow_num;
                overflow_key_num += ((InnerNode*)n)->Size();
            }
        }
        else if(dataslot[i].key!=FREE_FLAG && GET_TYPE(dataslot[i].val)==VALUE)
            ++key_num;
    }
    // key_num = this->key_num;
    fill_rate = (double)key_num / this->key_num;
    overflow_rate = (double)overflow_key_num / this->key_num;
}

uint64_t InnerNode::Memory() const{
    uint64_t res = sizeof(InnerNode);
    res += (capacity_ * sizeof(Record<_key_t, _value_t>));
    for (size_t i = 0; i < capacity_; i ++ ){
        if(GET_TYPE(dataslot[i].val) == OFNODE){
            InnerNode *n = (InnerNode *)GET_POINTER(dataslot[i].val);
            res += n->Memory();
        }
    }
    return res;
}

    
} // namespace lisk
