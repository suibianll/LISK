

#ifndef __LISK_ALEXNODE_H__
#define __LISK_ALEXNODE_H__
#include <cassert>
#include <cstring>
#include <atomic>
#include <immintrin.h>
#include <sched.h>
#include "node.h"
#include "util.h"

namespace lisk {
class ALEXNode : public Node{

public:
    ALEXNode(_key_t *keys, _value_t *values, uint64_t num, uint16_t d = 0);

    ALEXNode();
    ~ALEXNode(){delete idx;}

    bool Insert(_key_t key, _value_t &val);

    bool Find(_key_t key, _value_t &val);

    bool Update(_key_t key, _value_t &val);
    
    bool Remove(_key_t key);

    bool Scan(_key_t startKey, uint64_t len, _value_t* &result);

    void Stat(uint64_t &key_num, uint64_t &overflow_key_num, double & fill_rate, double &overflow_rate, uint64_t &overflow_num);

    uint64_t Size() const {return 0;}

    uint64_t Memory() const {return idx->model_size();sizeof(Node); sizeof(alex::Alex<_key_t, _value_t>);}

private:    
    uint32_t capacity_;
    
    alex::Alex<_key_t, _value_t>* idx;
	uint8_t buffer[226];
 
}__attribute((packed));
}//namespace lisk
#endif