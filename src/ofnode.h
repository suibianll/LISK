
#ifndef _LISK_OFNODE_H_
#define _LISK_OFNODE_H_

#include <string.h>
#include "util.h"
#include "hashnode.h"

namespace lisk{
#define SIMPLE_LEAF_NODE_LEVEL1_SIZE 12
#define SIMPLE_LEAF_NODE_LEVEL2_SIZE (SIMPLE_LEAF_NODE_LEVEL1_SIZE << 3)
#define SIMPLE_LEAF_NODE_LEVEL3_SIZE (SIMPLE_LEAF_NODE_LEVEL1_SIZE << 6)
#define Equal(k1, k2) (k1 == k2)
#define Greater_than(k1, k2) (k1 > k2)
#define Less_than(k1, k2) (k1 < k2)
extern uint64_t flat_probe_length;
extern uint64_t expand_failure;

template<class T, class P>
class OFArray{
public:
    uint32_t count_;
    uint32_t capacity_;
    Record<T, P> *dataslot_;

    OFArray(uint32_t cap){
        capacity_ = cap;
        count_ = 0;
        dataslot_ = new Record<T,P>[cap];
        // sizeof(HashMap<_key_t,_value_t>);
    }
    OFArray(T* keys, P* values, int num){
        capacity_ = num * 2;
        count_ = num;

        dataslot_ = new Record<T,P>[capacity_];
        for (size_t i = 0; i < num; i++){
            dataslot_[i] = {keys[i], values[i]};
        }
        
    }
    ~OFArray(){}

    inline bool is_Full(){return count_ == capacity_;}//for expanding the capacity
    
    inline bool IsFull(){return false;}//for outer calling

    inline bool IsEmpty(){return count_ == 0;};

    uint64_t Size(){
        return count_;
    }

    uint64_t Memory(){
        return sizeof(OFArray<T, P>) + capacity_ * (sizeof(T) + sizeof(P));
    }

    bool Find(const T &key, P &val){
        #ifndef NDEBUG
        of_search_time++;
        #endif
        if(count_ < 8){

            for (size_t i = 0; i < count_; i++){
                #ifndef NDEBUG
                of_probe_length++;
                #endif
                if(dataslot_[i].key == key){
                    val = dataslot_[i].val;
                    return true;
                }

            }
            return false;
        }
        else{
            
            int pos = BinarySearch(dataslot_, count_, key);

            if(pos != -1) {
                val = dataslot_[pos].val;
                return true;
            }
            return false;
        }
    }

    bool Locate(const T &key, P* &val){
        #ifndef NDEBUG
        of_search_time++;
        #endif
        if(count_ < 8){

            for (size_t i = 0; i < count_; i++){
                #ifndef NDEBUG
                of_probe_length++;
                #endif
                if(dataslot_[i].key == key){
                    val = &(dataslot_[i].val);
                    return true;
                }

            }
            return false;
        }
        else{
            
            int pos = BinarySearch(dataslot_, count_, key);

            if(pos != -1) {
                val = &(dataslot_[pos].val);
                return true;
            }
            return false;
        }
    }

    bool Insert(const T &key, const P &val){
        if(count_ == capacity_){//full
            capacity_ = 2 * capacity_;
            Record<T, P> *old_dataslot_ = dataslot_;
            Record<T, P> *new_dataslot_ = new Record<T,P>[capacity_];
            memcpy(new_dataslot_, old_dataslot_, count_ * sizeof(Record<T,P>));
            dataslot_ = new_dataslot_;
            delete[] old_dataslot_;
        }
        int pos = BinarySearch_LowerBound<T,P>(dataslot_, count_, key);
        memmove(&dataslot_[pos + 1], &dataslot_[pos], (count_ - pos) * sizeof(Record<T, P>));
        dataslot_[pos] = Record<T,P>(key, val);
        count_ ++;
        return true;
    }

    uint64_t FindInsert(const T &key, P &val){
        return 0;
    }

    bool Update(const T &key, const P &val){
        if(count_ < 8){
            for (size_t i = 0; i < count_; i++){
                if(dataslot_[i].key == key){
                    dataslot_[i].val = val;
                    return true;
                }
            }
            return false;
        }
        else{
            int pos = BinarySearch(dataslot_, count_, key);
            if(pos != -1){
                dataslot_[pos].val = val;
                return true;
            }
            return false;
        }
    }

    bool Remove(const T &key){
        if(count_ < 8){
            for (size_t i = 0; i < count_; i++){
                if(dataslot_[i].key == key){
                    memmove(&dataslot_[i], &dataslot_[i + 1], (count_ - 1 - i) * sizeof(Record<T,P>));
                    return true;
                }
            }
            return false;
        }
        else{
            int pos = BinarySearch(dataslot_, count_, key);
            if(pos != -1){
                memmove(&dataslot_[pos], &dataslot_[pos + 1], (count_ - 1 - pos) * sizeof(Record<T,P>));
                return true;
            }
            return false;
        }
    }

    int Scan(T &startKey, int len, Record<T,P> *result){
        
        int pos = BinarySearch_UpperBound(dataslot_, count_, startKey);
        int count = 0;
        for (size_t i = pos; i < count; i++){
            result[count++] = dataslot_[i];
            if(count == len) break;
        }
        return count;
        
    }
    
    void Dump(std::vector<_key_t> &keys, std::vector<_value_t> &values){
        for (size_t i = 0; i < count_; i++)
        {
            keys.emplace_back(dataslot_[i].key);
            values.emplace_back(dataslot_[i].val);
        }
        
    }

    void Dump(_key_t *keys, _value_t* values, uint64_t &start){
        for (size_t i = 0; i < count_; i++){
            keys[start] = dataslot_[i].key;
            values[start] = dataslot_[i].val;
            start ++;
        }
    }

    inline Record<T, P> PopFront(){
        Record<T, P> result = dataslot_[0];
        memmove(&dataslot_[0], &dataslot_[1], (count_ - 1) * sizeof(Record<T,P>));
        count_--;
        return result;
    }
    inline Record<T, P> Front(){
        return dataslot_[0];
    }
    __always_inline void Prefetch(){
        prefetch(this);
    }
};

// Left buffer, right buffer, overflow block.
template<typename T, typename P>
class FlatOFNode {


protected:

    uint8_t keys_level = 1;
    uint16_t number = 0;
    P* payloads;
    T keys_l1[SIMPLE_LEAF_NODE_LEVEL1_SIZE];
    T* keys_l2l3 = NULL;
static constexpr T free_flag = std::numeric_limits<T>::max();
    // Expand level 1 to level 2 (x8)
    void expand_l2() {
        keys_l2l3 = static_cast<T*>(std::aligned_alloc(64, sizeof(T) * SIMPLE_LEAF_NODE_LEVEL2_SIZE));
        P* payloads_tmp = payloads;
        payloads = new P[SIMPLE_LEAF_NODE_LEVEL2_SIZE];
        for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
            keys_l2l3[i << 3] = keys_l1[i];
            payloads[i << 3] = payloads_tmp[i];
            for (size_t j = 1; j < 8; j++) {
                keys_l2l3[(i << 3) + j] = free_flag;
            }
        }
        delete [] payloads_tmp;
        keys_level = 2;
    }

    // Expand level 2 to level 3 (x8)
    void expand_l3() {
        T* keys_l2_tmp = keys_l2l3;
        keys_l2l3 = static_cast<T*>(std::aligned_alloc(64, sizeof(T) * (SIMPLE_LEAF_NODE_LEVEL2_SIZE + SIMPLE_LEAF_NODE_LEVEL3_SIZE)));
        memcpy(keys_l2l3, keys_l2_tmp, SIMPLE_LEAF_NODE_LEVEL2_SIZE * sizeof(T));
        P* payloads_tmp = payloads;
        payloads = new P[SIMPLE_LEAF_NODE_LEVEL3_SIZE];
        for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++) {
            keys_l2l3[(i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = keys_l2_tmp[i];
            payloads[i << 3] = payloads_tmp[i];
            for (size_t j = 1; j < 8; j++) {
                keys_l2l3[(i << 3) + j + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = free_flag;
            }
        }
        free(keys_l2_tmp);
        delete [] payloads_tmp;
        keys_level = 3;
    }

public:
    // uint64_t count_;
    explicit FlatOFNode(int num):payloads(new P[SIMPLE_LEAF_NODE_LEVEL1_SIZE]) {
                for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                    keys_l1[i] = free_flag;
                }
            }
    explicit FlatOFNode(T* keys, P* values, int num){
        bulk_load(keys, values, num);        
    }
    ~FlatOFNode() {
        if (keys_level > 1) {
            assert(keys_l2l3 != NULL);
            free(keys_l2l3);
        }
        delete [] payloads;
    }

    __always_inline void Prefetch(){
        prefetch(&keys_level);
        prefetch(keys_l1);
    }

    void delete_children() { delete this; }

    size_t get_number() { return number; }

    T* get_keys() {
        if (keys_level == 1)
            return keys_l1;
        else if (keys_level == 2)
            return keys_l2l3;
        else
            return keys_l2l3 + SIMPLE_LEAF_NODE_LEVEL2_SIZE;
    }

    P* get_payloads() { return payloads; }

    size_t get_size() {
        if (keys_level == 1)
            return number;
        else if (keys_level == 2)
            return SIMPLE_LEAF_NODE_LEVEL2_SIZE;
        else
            return SIMPLE_LEAF_NODE_LEVEL3_SIZE;
    }

    bool bulk_load(T* keys, P* values, int num) {
        number = num;
        if(num < SIMPLE_LEAF_NODE_LEVEL1_SIZE){
            keys_level = 1;
            payloads = new P[SIMPLE_LEAF_NODE_LEVEL1_SIZE];
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                if(i < num){
                    keys_l1[i] = keys[i];
                    payloads[i] = values[i];
                }
                else{
                    keys_l1[i] = free_flag;
                }
                
            }
            return true;
        }
        else if(num < SIMPLE_LEAF_NODE_LEVEL2_SIZE){
            keys_level = 2;
            keys_l2l3 = static_cast<T*>(std::aligned_alloc(64, sizeof(T) * SIMPLE_LEAF_NODE_LEVEL2_SIZE));
            payloads = new P[SIMPLE_LEAF_NODE_LEVEL2_SIZE];
            int leaf_num = SIMPLE_LEAF_NODE_LEVEL1_SIZE;
            int leaf_size = num / leaf_num;
            int cur = 0;
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                keys_l1[i] = keys[cur];
                leaf_size = num / leaf_num;
                for (size_t j = 0; j < 8; j++) {
                    if(j < leaf_size){
                        keys_l2l3[(i << 3) + j] = keys[cur];
                        payloads[(i << 3) + j] = values[cur];
                        cur++;    
                    }
                    else{
                        keys_l2l3[(i << 3) + j] = free_flag;
                    }
                    
                }
                num-= leaf_size;
                leaf_num--;
            }
            return true;
        }
        else if(num <= SIMPLE_LEAF_NODE_LEVEL3_SIZE){
            keys_level = 3;
            keys_l2l3 = static_cast<T*>(std::aligned_alloc(64, sizeof(T) * (SIMPLE_LEAF_NODE_LEVEL2_SIZE + SIMPLE_LEAF_NODE_LEVEL3_SIZE)));
            payloads = new P[SIMPLE_LEAF_NODE_LEVEL3_SIZE];
            int leaf_num = SIMPLE_LEAF_NODE_LEVEL2_SIZE;
            int leaf_size = num / leaf_num;
            int cur = 0;
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++){
                keys_l2l3[i] = keys[cur];
                leaf_size = num / leaf_num;
                for (size_t j = 0; j < 8; j++) {
                    if(j < leaf_size){
                        keys_l2l3[(i << 3) + j + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = keys[cur];
                        payloads[(i << 3) + j] = values[cur];    
                        cur++;
                    }
                    else{
                        keys_l2l3[(i << 3) + j + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = free_flag;
                    }
                    
                }
                num-= leaf_size;
                leaf_num--;
            }
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                keys_l1[i] = keys_l2l3[i << 3];
            }
            
            return true;
        }
        else return false;
        
    }

    bool Find(const T &key, P &val) const {
        #ifndef NDEBUG
        of_search_time++;
        #endif
        // prefetch(keys_l1);

        // Level 1
        if (likely(keys_level == 1)) {
            prefetch(&payloads[0]);
            for (size_t i = 0; i < number; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Equal(keys_l1[i], key)){
                    val = payloads[i];
                    return true;
                }

            }
            return false;
        }

        // Level 2
        else if (keys_level == 2) {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            prefetch(&payloads[start_pos_l2]);
            if (start_pos_l2 >= 0){
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l2 + i], key)){
                        val = payloads[start_pos_l2 + i];
                        return true;
                    }
                }
            }
            return false;
        }

        // Level 3
        else {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0) {
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (size_t i = 1; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Greater_than(keys_l2l3[start_pos_l2 + i], key)) {
                        start_pos_l3 = (start_pos_l2 + i - 1) << 3;
                        break;
                    }
                }
                prefetch(&payloads[start_pos_l3]);
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)){
                        val = payloads[start_pos_l3 + i];
                        return true;
                    }
                }
            }
            return false;
        }
    }

    bool Locate(const T &key, P* &val) const {
        #ifndef NDEBUG
        of_search_time++;
        #endif

        // Level 1
        if (keys_level == 1) {
            prefetch(&payloads[0]);
            for (size_t i = 0; i < number; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Equal(keys_l1[i], key)){
                    val = &payloads[i];
                    return true;
                }

            }
            return false;
        }

        // Level 2
        else if (keys_level == 2) {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0){
                prefetch(&payloads[start_pos_l2]);
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l2 + i], key)){
                        val = &payloads[start_pos_l2 + i];
                        return true;
                    }
                }
            }
            return false;
        }

        // Level 3
        else {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0) {
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (size_t i = 1; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Greater_than(keys_l2l3[start_pos_l2 + i], key)) {
                        start_pos_l3 = (start_pos_l2 + i - 1) << 3;
                        break;
                    }
                }
                prefetch(&payloads[start_pos_l3]);
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)){
                        val = &payloads[start_pos_l3 + i];
                        return true;
                    }
                }
            }
            return false;
        }
    }

    void range_query(bool type, T lower_bound, T upper_bound, std::vector<std::pair<T, P>>& answers) const {
        if (type) {
            // Level 1
            if (keys_level == 1) {
                for (size_t it = 0; it < number; it++) {
                    if (!Greater_than(keys_l1[it], upper_bound))
                        answers.emplace_back(keys_l1[it], payloads[it]);
                    else
                        return;
                }
            }

            // Level 2
            else if (keys_level == 2) {
                for (size_t it = 0; it < SIMPLE_LEAF_NODE_LEVEL2_SIZE; it++) {
                    if(!Greater_than(keys_l2l3[it], upper_bound)) {
                        if (!Equal(keys_l2l3[it], free_flag))
                            answers.emplace_back(keys_l2l3[it], payloads[it]);
                    }
                    else
                        return;
                }
            }

            // Level 3
            else {
                for (size_t it = 0; it < SIMPLE_LEAF_NODE_LEVEL3_SIZE; it++) {
                    if (!Greater_than(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], upper_bound)) {
                        if (!Equal(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag))
                            answers.emplace_back(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], payloads[it]);
                    }
                    else
                        return;
                }
            }
        }
        else {
            bool end_flag = false;
            // Level 1
            if (keys_level == 1) {
                for (size_t it = 0; it < number; it++) {
                    if (!Less_than(keys_l1[it], lower_bound)) {
                        while (it < number && !Greater_than(keys_l1[it], upper_bound)) {
                            answers.emplace_back(keys_l1[it], payloads[it]);
                            it++;
                        }
                        if (it < number)
                            end_flag = true;
                        break;
                    }
                }
            }

            // Level 2
            else if (keys_level == 2) {
                int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
                for (size_t i = 1; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                    if (Greater_than(keys_l1[i], lower_bound)) {
                        start_pos_l2 = (i - 1) << 3;
                        break;
                    }
                }
                for (size_t i = 0; i < 8; i++) {
                    if (!Less_than(keys_l2l3[start_pos_l2 + i], lower_bound)) {
                        size_t it = start_pos_l2 + i;
                        while (it < SIMPLE_LEAF_NODE_LEVEL2_SIZE && !Greater_than(keys_l2l3[it], upper_bound)) {
                            if (!Equal(keys_l2l3[it], free_flag))
                                answers.emplace_back(keys_l2l3[it], payloads[it]);
                            it++;
                        }
                        if (it < SIMPLE_LEAF_NODE_LEVEL2_SIZE)
                            end_flag = true;
                        break;
                    }
                }
            }

            // Level 3
            else {
                int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
                for (size_t i = 1; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                    if (Greater_than(keys_l1[i], lower_bound)) {
                        start_pos_l2 = (i - 1) << 3;
                        break;
                    }
                }
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (size_t i = 1; i < 8; i++) {
                    if (Greater_than(keys_l2l3[start_pos_l2 + i], lower_bound)) {
                        start_pos_l3 = (start_pos_l2 + i - 1) << 3;
                        break;
                    }
                }
                for (size_t i = 0; i < 8; i++) {
                    if (!Less_than(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], lower_bound)) {
                        size_t it = start_pos_l3 + i;
                        while (it < SIMPLE_LEAF_NODE_LEVEL3_SIZE && !Greater_than(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], upper_bound)) {
                            if (!Equal(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag))
                                answers.emplace_back(keys_l2l3[it + SIMPLE_LEAF_NODE_LEVEL2_SIZE], payloads[it]);
                            it++;
                        }
                        if (it < SIMPLE_LEAF_NODE_LEVEL3_SIZE)
                            end_flag = true;
                        break;
                    }
                }
            }

        }
    }

    // void find_split_node(T key, std::vector<NODE*>& traversal_path) {
    //     throw std::logic_error("Get traversal path error!");
    // }

    // 1: insert success; 2: update success; -1: insert fault and evolve node.
    bool Insert(const T &key, const P &val) {
        // assert(number < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
        // Level 1
        if (keys_level == 1) {
            size_t target_pos = number;
            for (size_t i = 0; i < number; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    target_pos = i;
                    break;
                }
            }
            // Update
            if (target_pos > 0 && Equal(keys_l1[target_pos - 1], key)) {
                payloads[target_pos - 1] = val;
                return true;
            }
            assert(number >= target_pos && number < SIMPLE_LEAF_NODE_LEVEL1_SIZE);
            memmove(keys_l1 + target_pos + 1, keys_l1 + target_pos, (number - target_pos) * sizeof(T));
            memmove(payloads + target_pos + 1, payloads + target_pos, (number - target_pos) * sizeof(P));
            keys_l1[target_pos] = key;
            payloads[target_pos] = val;
            number++;
            if (number == SIMPLE_LEAF_NODE_LEVEL1_SIZE)
                expand_l2();
            return true;
        }

        // Level 2
        else if (keys_level == 2) {

            // Get start position (line) at level 2
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0) {
                int target_pos = start_pos_l2 + 8;
                int free_pos = target_pos;
                for (size_t i = 0; i < 8; i++) {
                    // Update
                    if (Equal(keys_l2l3[start_pos_l2 + i], key)) {
                        payloads[start_pos_l2 + i] = val;
                        return true;
                    }
                    else if (Equal(keys_l2l3[start_pos_l2 + i], free_flag)) {
                        free_pos = start_pos_l2 + i;
                        target_pos = free_pos;
                        break;
                    }
                    else if (Greater_than(keys_l2l3[start_pos_l2 + i], key)) {
                        target_pos = start_pos_l2 + i;
                        break;
                    }
                }
                if (free_pos != target_pos) {
                    for (int i = target_pos - start_pos_l2 + 1; i < 8; i++) {
                        if (Equal(keys_l2l3[start_pos_l2 + i], free_flag)) {
                            free_pos = start_pos_l2 + i;
                            break;
                        }
                    }
                }
                // Target line has free position
                if (free_pos < start_pos_l2 + 8) {
                    assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                    memmove(keys_l2l3 + target_pos + 1, keys_l2l3 + target_pos, (free_pos - target_pos) * sizeof(T));
                    memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                    keys_l2l3[target_pos] = key;
                    payloads[target_pos] = val;
                }
                // Target line is full
                else {
                    for (int i = 1; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                        // Borrow from front line
                        if (start_pos_l2 - (i << 3) + 7 > 0 && Equal(keys_l2l3[start_pos_l2 - (i << 3) + 7], free_flag)) {
                            size_t influenced_header = start_pos_l2 - (i << 3) + 8;
                            for (size_t j = 2; j < 9; j++) {
                                if (!Equal(keys_l2l3[influenced_header - j], free_flag)) {
                                    free_pos = influenced_header - j + 1;
                                    break;
                                }
                            }
                            keys_l2l3[free_pos] = keys_l2l3[influenced_header];
                            payloads[free_pos] = payloads[influenced_header];
                            assert(target_pos > influenced_header && target_pos <= SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                            memmove(keys_l2l3 + influenced_header, keys_l2l3 + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(T));
                            memmove(payloads + influenced_header, payloads + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(P));
                            keys_l2l3[target_pos - 1] = key;
                            payloads[target_pos - 1] = val;
                            for (size_t j = influenced_header; j < target_pos; j = j + 8) {
                                keys_l1[j >> 3] = keys_l2l3[j];
                            }
                            break;
                        }
                        // Borrow from latter line
                        if (start_pos_l2 + (i << 3) + 7 < SIMPLE_LEAF_NODE_LEVEL2_SIZE && Equal(keys_l2l3[start_pos_l2 + (i << 3) + 7], free_flag)) {
                            size_t influenced_header = start_pos_l2 + (i << 3);
                            for (size_t j = 1; j < 8; j++) {
                                if (Equal(keys_l2l3[influenced_header + j], free_flag)) {
                                    free_pos = influenced_header + j;
                                    break;
                                }
                            }
                            assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                            memmove(keys_l2l3 + target_pos + 1, keys_l2l3 + target_pos, (free_pos - target_pos) * sizeof(T));
                            memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                            keys_l2l3[target_pos] = key;
                            payloads[target_pos] = val;
                            for (int j = influenced_header; j >= target_pos; j = j - 8) {
                                keys_l1[j >> 3] = keys_l2l3[j];
                            }
                            break;
                        }
                    }
                }
            }
            // Insert in the first position
            else {
                size_t free_pos = 7;
                for (size_t i = 7; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i = i + 8) {
                    if (Equal(keys_l2l3[i], free_flag)) {
                        free_pos = i;
                        break;
                    }
                }
                for (size_t i = 1; i < 8; i++) {
                    if (!Equal(keys_l2l3[free_pos - i], free_flag)) {
                        free_pos = free_pos - i + 1;
                        break;
                    }
                }
                assert(free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                memmove(keys_l2l3 + 1, keys_l2l3, free_pos * sizeof(T));
                memmove(payloads + 1, payloads, free_pos * sizeof(P));
                keys_l2l3[0] = key;
                payloads[0] = val;
                for (size_t i = 0; i << 3 < free_pos; i++) {
                    keys_l1[i] = keys_l2l3[i << 3];
                }
            }
            number++;
            if (number == SIMPLE_LEAF_NODE_LEVEL2_SIZE) {
                expand_l3();
            }
                
            return true;
        }

        // Level 3
        else {
            // Get start position (line) at level 2
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }

            if (start_pos_l2 >= 0) {
                
                // Get start position (line) at level 3
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (int i = start_pos_l2; i < start_pos_l2 + 8; i++) {
                    if (Greater_than(keys_l2l3[i], key)) {
                        start_pos_l3 = (i - 1) << 3;
                        break;
                    }
                }

                int target_pos = start_pos_l3 + 8;
                int free_pos = target_pos;
                for (size_t i = 0; i < 8; i++) {
                    // Update
                    if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)) {
                        payloads[start_pos_l3 + i] = val;
                        return true;
                    }
                    else if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = start_pos_l3 + i;
                        target_pos = free_pos;
                        break;
                    }
                    else if (Greater_than(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)) {
                        target_pos = start_pos_l3 + i;
                        break;
                    }
                }
                if (free_pos != target_pos) {
                    for (int i = target_pos - start_pos_l3 + 1; i < 8; i++) {
                        if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                            free_pos = start_pos_l3 + i;
                            break;
                        }
                    }
                }
                // Target line has free position
                if (free_pos < start_pos_l3 + 8) {
                    assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                    memmove(keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE, (free_pos - target_pos) * sizeof(T));
                    memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                    keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                    payloads[target_pos] = val;
                }
                // Borrow free position
                else {
                    for (int i = 1; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++) {
                        // Borrow from front line
                        if (start_pos_l3 - (i << 3) + 7 > 0 && Equal(keys_l2l3[start_pos_l3 - (i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 7], free_flag)) {
                            size_t influenced_header = start_pos_l3 - (i << 3) + 8;
                            for (size_t j = 2; j < 9; j++) {
                                if (!Equal(keys_l2l3[influenced_header - j + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                                    free_pos = influenced_header - j + 1;
                                    break;
                                }
                            }
                            keys_l2l3[free_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = keys_l2l3[influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            payloads[free_pos] = payloads[influenced_header];
                            assert(target_pos > influenced_header && target_pos <= SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                            memmove(keys_l2l3 + influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE, keys_l2l3 + influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, (target_pos - influenced_header - 1) * sizeof(T));
                            memmove(payloads + influenced_header, payloads + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(P));
                            keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE - 1] = key;
                            payloads[target_pos - 1] = val;
                            for (size_t j = influenced_header; j < target_pos; j = j + 8) {
                                keys_l2l3[j >> 3] = keys_l2l3[j + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            }
                            for (size_t j = ((influenced_header - 1) >> 6) + 1; j << 6 < target_pos; j++) {
                                keys_l1[j] = keys_l2l3[j << 3];
                            }
                            break;
                        }
                        // Borrow from latter line
                        if (start_pos_l3 + (i << 3) + 7 < SIMPLE_LEAF_NODE_LEVEL3_SIZE && Equal(keys_l2l3[start_pos_l3 + (i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 7], free_flag)) {
                            size_t influenced_header = start_pos_l3 + (i << 3);
                            for (size_t j = 1; j < 8; j++) {
                                if (Equal(keys_l2l3[influenced_header + j + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                                    free_pos = influenced_header + j;
                                    break;
                                }
                            }
                            assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                            memmove(keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE, (free_pos - target_pos) * sizeof(T));
                            memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                            keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                            payloads[target_pos] = val;
                            for (int j = influenced_header; j >= target_pos; j = j - 8) {
                                keys_l2l3[j >> 3] = keys_l2l3[j + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            }
                            for (int j = influenced_header >> 6; j << 6 >= target_pos; j--) {
                                keys_l1[j] = keys_l2l3[j << 3];
                            }
                            break;
                        }
                    }
                }
            }
            // Insert in the first position
            else {
                size_t free_pos = 7;
                for (size_t i = 7; i < SIMPLE_LEAF_NODE_LEVEL3_SIZE; i = i + 8) {
                    if (Equal(keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = i;
                        break;
                    }
                }
                for (size_t i = 1; i < 8; i++) {
                    if (!Equal(keys_l2l3[free_pos - i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = free_pos - i + 1;
                        break;
                    }
                }
                assert(free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                memmove(keys_l2l3 + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + SIMPLE_LEAF_NODE_LEVEL2_SIZE, free_pos * sizeof(T));
                memmove(payloads + 1, payloads, free_pos * sizeof(P));
                keys_l2l3[SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                payloads[0] = val;
                for (size_t i = 0; i << 3 < free_pos; i++) {
                    keys_l2l3[i] = keys_l2l3[(i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                }
                for (size_t i = 0; i << 6 < free_pos; i++) {
                    keys_l1[i] = keys_l2l3[i << 3];
                }
            }

            number++;
            // Evolve into learned node
            if (number <= SIMPLE_LEAF_NODE_LEVEL2_SIZE * 6)
                return true;
            else
                return false;
        }
    }

    uint64_t FindInsert(const T &key, const P &val) {//try to insert, if the key already exists, return the address
        // assert(number < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
        // Level 1
        // LOG(INFO)<<"FindInsert: "<<key<<" oftree"<<this<<" key num:"<<number<<" payloads address"<<payloads;
        if (keys_level == 1) {
            size_t target_pos = number;
            for (size_t i = 0; i < number; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    target_pos = i;
                    break;
                }
            }
            // Update
            if (target_pos > 0 && Equal(keys_l1[target_pos - 1], key)) {
                // payloads[target_pos - 1] = val;
                return uint64_t(&payloads[target_pos - 1]);
            }
            assert(number >= target_pos && number < SIMPLE_LEAF_NODE_LEVEL1_SIZE);
            memmove(keys_l1 + target_pos + 1, keys_l1 + target_pos, (number - target_pos) * sizeof(T));
            memmove(payloads + target_pos + 1, payloads + target_pos, (number - target_pos) * sizeof(P));
            keys_l1[target_pos] = key;
            payloads[target_pos] = val;
            number++;
            if (number == SIMPLE_LEAF_NODE_LEVEL1_SIZE)
                expand_l2();
            return 0;
        }

        // Level 2
        else if (keys_level == 2) {

            // Get start position (line) at level 2
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0) {
                int target_pos = start_pos_l2 + 8;
                int free_pos = target_pos;
                for (size_t i = 0; i < 8; i++) {
                    // Update
                    if (Equal(keys_l2l3[start_pos_l2 + i], key)) {
                        // payloads[start_pos_l2 + i] = val;
                        return uint64_t(&payloads[start_pos_l2 + i]);
                    }
                    else if (Equal(keys_l2l3[start_pos_l2 + i], free_flag)) {
                        free_pos = start_pos_l2 + i;
                        target_pos = free_pos;
                        break;
                    }
                    else if (Greater_than(keys_l2l3[start_pos_l2 + i], key)) {
                        target_pos = start_pos_l2 + i;
                        break;
                    }
                }
                if (free_pos != target_pos) {
                    for (int i = target_pos - start_pos_l2 + 1; i < 8; i++) {
                        if (Equal(keys_l2l3[start_pos_l2 + i], free_flag)) {
                            free_pos = start_pos_l2 + i;
                            break;
                        }
                    }
                }
                // Target line has free position
                if (free_pos < start_pos_l2 + 8) {
                    assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                    memmove(keys_l2l3 + target_pos + 1, keys_l2l3 + target_pos, (free_pos - target_pos) * sizeof(T));
                    memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                    keys_l2l3[target_pos] = key;
                    payloads[target_pos] = val;
                }
                // Target line is full
                else {
                    for (int i = 1; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                        // Borrow from front line
                        if (start_pos_l2 - (i << 3) + 7 > 0 && Equal(keys_l2l3[start_pos_l2 - (i << 3) + 7], free_flag)) {
                            size_t influenced_header = start_pos_l2 - (i << 3) + 8;
                            for (size_t j = 2; j < 9; j++) {
                                if (!Equal(keys_l2l3[influenced_header - j], free_flag)) {
                                    free_pos = influenced_header - j + 1;
                                    break;
                                }
                            }
                            keys_l2l3[free_pos] = keys_l2l3[influenced_header];
                            payloads[free_pos] = payloads[influenced_header];
                            assert(target_pos > influenced_header && target_pos <= SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                            memmove(keys_l2l3 + influenced_header, keys_l2l3 + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(T));
                            memmove(payloads + influenced_header, payloads + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(P));
                            keys_l2l3[target_pos - 1] = key;
                            payloads[target_pos - 1] = val;
                            for (size_t j = influenced_header; j < target_pos; j = j + 8) {
                                keys_l1[j >> 3] = keys_l2l3[j];
                            }
                            break;
                        }
                        // Borrow from latter line
                        if (start_pos_l2 + (i << 3) + 7 < SIMPLE_LEAF_NODE_LEVEL2_SIZE && Equal(keys_l2l3[start_pos_l2 + (i << 3) + 7], free_flag)) {
                            size_t influenced_header = start_pos_l2 + (i << 3);
                            for (size_t j = 1; j < 8; j++) {
                                if (Equal(keys_l2l3[influenced_header + j], free_flag)) {
                                    free_pos = influenced_header + j;
                                    break;
                                }
                            }
                            assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                            memmove(keys_l2l3 + target_pos + 1, keys_l2l3 + target_pos, (free_pos - target_pos) * sizeof(T));
                            memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                            keys_l2l3[target_pos] = key;
                            payloads[target_pos] = val;
                            for (int j = influenced_header; j >= target_pos; j = j - 8) {
                                keys_l1[j >> 3] = keys_l2l3[j];
                            }
                            break;
                        }
                    }
                }
            }
            // Insert in the first position
            else {
                size_t free_pos = 7;
                for (size_t i = 7; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i = i + 8) {
                    if (Equal(keys_l2l3[i], free_flag)) {
                        free_pos = i;
                        break;
                    }
                }
                for (size_t i = 1; i < 8; i++) {
                    if (!Equal(keys_l2l3[free_pos - i], free_flag)) {
                        free_pos = free_pos - i + 1;
                        break;
                    }
                }
                assert(free_pos < SIMPLE_LEAF_NODE_LEVEL2_SIZE);
                memmove(keys_l2l3 + 1, keys_l2l3, free_pos * sizeof(T));
                memmove(payloads + 1, payloads, free_pos * sizeof(P));
                keys_l2l3[0] = key;
                payloads[0] = val;
                for (size_t i = 0; i << 3 < free_pos; i++) {
                    keys_l1[i] = keys_l2l3[i << 3];
                }
            }
            number++;
            if (number == SIMPLE_LEAF_NODE_LEVEL2_SIZE) {
                expand_l3();
            }
                
            return 0;
        }

        // Level 3
        else {
            // Get start position (line) at level 2
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }

            if (start_pos_l2 >= 0) {
                
                // Get start position (line) at level 3
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (int i = start_pos_l2; i < start_pos_l2 + 8; i++) {
                    if (Greater_than(keys_l2l3[i], key)) {
                        start_pos_l3 = (i - 1) << 3;
                        break;
                    }
                }

                int target_pos = start_pos_l3 + 8;
                int free_pos = target_pos;
                for (size_t i = 0; i < 8; i++) {
                    // Update
                    if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)) {
                        // payloads[start_pos_l3 + i] = val;
                        return uint64_t(&payloads[start_pos_l3 + i]);
                    }
                    else if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = start_pos_l3 + i;
                        target_pos = free_pos;
                        break;
                    }
                    else if (Greater_than(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)) {
                        target_pos = start_pos_l3 + i;
                        break;
                    }
                }
                if (free_pos != target_pos) {
                    for (int i = target_pos - start_pos_l3 + 1; i < 8; i++) {
                        if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                            free_pos = start_pos_l3 + i;
                            break;
                        }
                    }
                }
                // Target line has free position
                if (free_pos < start_pos_l3 + 8) {
                    assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                    memmove(keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE, (free_pos - target_pos) * sizeof(T));
                    memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                    keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                    payloads[target_pos] = val;
                }
                // Borrow free position
                else {
                    for (int i = 1; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++) {
                        // Borrow from front line
                        if (start_pos_l3 - (i << 3) + 7 > 0 && Equal(keys_l2l3[start_pos_l3 - (i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 7], free_flag)) {
                            size_t influenced_header = start_pos_l3 - (i << 3) + 8;
                            for (size_t j = 2; j < 9; j++) {
                                if (!Equal(keys_l2l3[influenced_header - j + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                                    free_pos = influenced_header - j + 1;
                                    break;
                                }
                            }
                            keys_l2l3[free_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = keys_l2l3[influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            payloads[free_pos] = payloads[influenced_header];
                            assert(target_pos > influenced_header && target_pos <= SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                            memmove(keys_l2l3 + influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE, keys_l2l3 + influenced_header + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, (target_pos - influenced_header - 1) * sizeof(T));
                            memmove(payloads + influenced_header, payloads + influenced_header + 1, (target_pos - influenced_header - 1) * sizeof(P));
                            keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE - 1] = key;
                            payloads[target_pos - 1] = val;
                            for (size_t j = influenced_header; j < target_pos; j = j + 8) {
                                keys_l2l3[j >> 3] = keys_l2l3[j + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            }
                            for (size_t j = ((influenced_header - 1) >> 6) + 1; j << 6 < target_pos; j++) {
                                keys_l1[j] = keys_l2l3[j << 3];
                            }
                            break;
                        }
                        // Borrow from latter line
                        if (start_pos_l3 + (i << 3) + 7 < SIMPLE_LEAF_NODE_LEVEL3_SIZE && Equal(keys_l2l3[start_pos_l3 + (i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 7], free_flag)) {
                            size_t influenced_header = start_pos_l3 + (i << 3);
                            for (size_t j = 1; j < 8; j++) {
                                if (Equal(keys_l2l3[influenced_header + j + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                                    free_pos = influenced_header + j;
                                    break;
                                }
                            }
                            assert(free_pos >= target_pos && free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                            memmove(keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE, (free_pos - target_pos) * sizeof(T));
                            memmove(payloads + target_pos + 1, payloads + target_pos, (free_pos - target_pos) * sizeof(P));
                            keys_l2l3[target_pos + SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                            payloads[target_pos] = val;
                            for (int j = influenced_header; j >= target_pos; j = j - 8) {
                                keys_l2l3[j >> 3] = keys_l2l3[j + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                            }
                            for (int j = influenced_header >> 6; j << 6 >= target_pos; j--) {
                                keys_l1[j] = keys_l2l3[j << 3];
                            }
                            break;
                        }
                    }
                }
            }
            // Insert in the first position
            else {
                size_t free_pos = 7;
                for (size_t i = 7; i < SIMPLE_LEAF_NODE_LEVEL3_SIZE; i = i + 8) {
                    if (Equal(keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = i;
                        break;
                    }
                }
                for (size_t i = 1; i < 8; i++) {
                    if (!Equal(keys_l2l3[free_pos - i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)) {
                        free_pos = free_pos - i + 1;
                        break;
                    }
                }
                assert(free_pos < SIMPLE_LEAF_NODE_LEVEL3_SIZE);
                memmove(keys_l2l3 + SIMPLE_LEAF_NODE_LEVEL2_SIZE + 1, keys_l2l3 + SIMPLE_LEAF_NODE_LEVEL2_SIZE, free_pos * sizeof(T));
                memmove(payloads + 1, payloads, free_pos * sizeof(P));
                keys_l2l3[SIMPLE_LEAF_NODE_LEVEL2_SIZE] = key;
                payloads[0] = val;
                for (size_t i = 0; i << 3 < free_pos; i++) {
                    keys_l2l3[i] = keys_l2l3[(i << 3) + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                }
                for (size_t i = 0; i << 6 < free_pos; i++) {
                    keys_l1[i] = keys_l2l3[i << 3];
                }
            }

            number++;
            return 0;
            // Evolve into learned node
            // if (number <= SIMPLE_LEAF_NODE_LEVEL2_SIZE * 6)
            //     return true;
            // else
            //     return false;
        }
    }

    void node_size(size_t& model_size, size_t& slot_size, size_t& data_size, size_t& inner_node, size_t& data_node) const {
        {
            model_size += sizeof(*this);
            data_size += number * (sizeof(T) + sizeof(P));
            if (keys_level == 1) {
                model_size -= SIMPLE_LEAF_NODE_LEVEL1_SIZE * sizeof(T);
                slot_size += SIMPLE_LEAF_NODE_LEVEL1_SIZE * (sizeof(T) + sizeof(P));
            }
            else if (keys_level == 2) {
                slot_size += SIMPLE_LEAF_NODE_LEVEL2_SIZE * (sizeof(T) + sizeof(P));
            }
            else {
                model_size += SIMPLE_LEAF_NODE_LEVEL2_SIZE * sizeof(T);
                slot_size += SIMPLE_LEAF_NODE_LEVEL3_SIZE * (sizeof(T) + sizeof(P));
            }
        }
    }

    bool Update(const T &key, const P &val){
        // Level 1
        if (likely(keys_level == 1)) {
            prefetch(&payloads[0]);
            for (size_t i = 0; i < number; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Equal(keys_l1[i], key)){
                    payloads[i] = val;
                    return true;
                }

            }
        }

        // Level 2
        else if (keys_level == 2) {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            prefetch(&payloads[start_pos_l2]);
            if (start_pos_l2 >= 0){
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l2 + i], key)){
                        payloads[start_pos_l2 + i] = val;
                        return true;
                    }
                }
            }
        }

        // Level 3
        else {
            int start_pos_l2 = SIMPLE_LEAF_NODE_LEVEL2_SIZE - 8;
            for (int i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++) {
                #ifndef NDEBUG
                    of_probe_length++;
                #endif
                if (Greater_than(keys_l1[i], key)) {
                    start_pos_l2 = (i - 1) << 3;
                    break;
                }
            }
            if (start_pos_l2 >= 0) {
                int start_pos_l3 = (start_pos_l2 << 3) + 56;
                for (size_t i = 1; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Greater_than(keys_l2l3[start_pos_l2 + i], key)) {
                        start_pos_l3 = (start_pos_l2 + i - 1) << 3;
                        break;
                    }
                }
                prefetch(&payloads[start_pos_l3]);
                for (size_t i = 0; i < 8; i++) {
                    #ifndef NDEBUG
                        of_probe_length++;
                    #endif
                    if (Equal(keys_l2l3[start_pos_l3 + i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], key)){
                        payloads[start_pos_l3 + i]  = val;
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool Remove(const T &key){
        return true;
    }

    int Scan(T &startKey, int len, Record<T,P> *result){
        return 0;
    }

    bool IsFull(){
        assert(number <= SIMPLE_LEAF_NODE_LEVEL3_SIZE);
        return number == SIMPLE_LEAF_NODE_LEVEL3_SIZE;
    }

    bool IsEmpty(){
        return number == 0;
    }
    
    void Dump(_key_t *keys, _value_t* values, uint64_t &start){//dump all the records in order to the vector
        if(keys_level == 1){
            for(size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                if(!Equal(keys_l1[i], free_flag)){
                    keys[start] = keys_l1[i];
                    values[start] = payloads[i];
                    start++;
                }
            }        
        }
        else if(keys_level == 2){
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++){
                if(!Equal(keys_l2l3[i], free_flag)){
                    keys[start] = keys_l2l3[i];
                    values[start] = payloads[i];
                    start ++;
                }
            }
            
        }
        else{
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL3_SIZE; i++){
                if(!Equal(keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag)){
                    keys[start] = keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE];
                    values[start] = payloads[i];
                    start ++;
                }
            }
        } 
    }

    void Print(){
        if(keys_level == 1){
            for(size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL1_SIZE; i++){
                // if(!Equal(keys_l1[i], free_flag))
                {
                    LOG(INFO)<<"flat node:"<<keys_l1[i]<<","<<payloads[i]<<"\n";
                }
            }        
        }
        else if(keys_level == 2){
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL2_SIZE; i++){
                // if(!Equal(keys_l2l3[i], free_flag))
                {
                    LOG(INFO)<<"flat node:"<<keys_l2l3[i]<<","<<payloads[i]<<"\n";
                }
            }
            
        }
        else{
            for (size_t i = 0; i < SIMPLE_LEAF_NODE_LEVEL3_SIZE; i++){
                // if(!Equal(keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], free_flag))
                {
                    LOG(INFO)<<"flat node:"<<keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE]<<","<<payloads[i]<<"\n";
                    // out.push_back(std::make_pair(keys_l2l3[i + SIMPLE_LEAF_NODE_LEVEL2_SIZE], payloads[i]));
                }
            }
        }
    }
    inline Record<T, P> PopFront(){
        return Record<T,P>();
    }
    inline Record<T, P> Front(){
        return Record<T,P>();
    }
    uint64_t Size(){
        return number;
    }

    uint64_t Memory(){
        uint64_t res = 0;
        res += sizeof(*this);
        if (keys_level == 1) {
            res += SIMPLE_LEAF_NODE_LEVEL1_SIZE * (sizeof(T) + sizeof(P));
        }
        else if (keys_level == 2) {
            res += SIMPLE_LEAF_NODE_LEVEL2_SIZE * (sizeof(T) + sizeof(P));
        }
        else {
            res += SIMPLE_LEAF_NODE_LEVEL2_SIZE * (sizeof(T));
            res += SIMPLE_LEAF_NODE_LEVEL3_SIZE * (sizeof(T) + sizeof(P));
        }
        return res;
    }
};//__attribute((packed));
}

#endif