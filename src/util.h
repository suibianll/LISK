
#ifndef __LISK_UTIL_H__
#define __LISK_UTIL_H__
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <assert.h>
#include <immintrin.h>
#include <glog/logging.h>
#include <queue>
#include "Key.h"

//statistics
// #define NDEBUG
#define NSTAT

#define SEP_ARRAYS
//define the type of the overflow node
// #define CLHT
#define HASHMAP
// #define ARRAY
// #define BTREE
// #define FLAT


namespace lisk{

const uint64_t SCALE_FACTOR = 2;
const _key_t MAX_KEY = std::numeric_limits<_key_t>::max();
const _key_t MIN_KEY = typeid(_key_t) == typeid(double) || typeid(_key_t) == typeid(float) 
                            ? -1 * MAX_KEY : std::numeric_limits<_key_t>::min();
const _key_t FREE_FLAG = std::numeric_limits<_key_t>::max();

extern uint64_t index_depth;
extern uint64_t inner_probe_length;
extern uint64_t inner_insert;
extern uint64_t leaf_probe_length;
extern uint64_t of_probe_length;
extern uint64_t of_search_time;
extern uint64_t retrain_time;
extern uint64_t hot_search;
extern uint64_t btree_search;
extern uint64_t learned_search;
extern uint64_t verify_num;
extern uint64_t total_leaf_num;
extern uint64_t fmcd_D;
extern uint64_t fmcd_success;
extern uint64_t fmcd_failure;
extern uint64_t expand_failure;
extern uint64_t search_depth[1000];
extern double latency_breakdown[100];
extern uint64_t bytes[8][256];

// #define likely(x)   __builtin_expect(!!(x), 1)
// #define unlikely(x) __builtin_expect(!!(x), 0)
#define CACHE_LINE_SIZE 64
typedef struct { char x[CACHE_LINE_SIZE]; } cacheline_t;
__always_inline void prefetch(const void *ptr) {
#ifdef NOPREFETCH
    (void) ptr;
#else
    asm volatile("prefetcht0 %0" : : "m" (*(const cacheline_t *)ptr));
#endif
}
enum TimeBreakdown{
  INNERTIME = 0,
  LEAFTIME,
  INNERPROBE,
  LEAFPROBE,
  OVERFLOWPROBE
};

const uint64_t POINTER_MASK = 0x0000ffffffffffff;
const uint64_t RADIX_MASK = 0x000fffffffffffff;
const uint64_t LOCK_MASK = 0x1000000000000000;
const uint64_t POINTER_TYPE_MASK = 0x0010000000000000;
const uint64_t POINTER_MASK_MOVE_NUM = 52;

enum NodeType {MINNER = 1, MLEAF, BINNER, BLEAF, SLEAF, SUBHOT, SUBART, SUBALEX};
enum PointerType {VALUE = 0, OFNODE = 5, DUPVALUE, DUPOFNODE, NEXTLEVEL};

#define SET_TYPE(x, type) ((uint64_t)x | ((uint64_t)type << POINTER_MASK_MOVE_NUM))
#define GET_TYPE(x) (x >> POINTER_MASK_MOVE_NUM)
#define GET_POINTER(x) ((uint64_t)x & POINTER_MASK)

template<class T, class P>
struct Record {
    T key;
    P val;
    
    Record(): key(std::numeric_limits<T>::max()), val(0) {}
    Record(T k, P v) : key(k), val(v) {}

    inline bool operator < (const Record & oth) {
        return key < oth.key;
    }

    inline bool operator > (const Record & oth) const {
        return key > oth.key;
    }

    inline void operator = (const Record & oth) {
        key = oth.key;
        val = oth.val;
    }
};//__attribute((packed));

template <class T>
class LinearModelBuilder {
public:
    double a_;
    double b_;
    LinearModelBuilder() { }

    inline void add(T x, int y) {
        count_++;
        x_sum_ += static_cast<long double>(x);
        y_sum_ += static_cast<long double>(y);
        xx_sum_ += static_cast<long double>(x) * x;
        xy_sum_ += static_cast<long double>(x) * y;
        x_min_ = std::min(x, x_min_);
        x_max_ = std::max(x, x_max_);
        y_min_ = std::min<double>(y, y_min_);
        y_max_ = std::max<double>(y, y_max_);
    }

    void build() {
        if (count_ <= 1) {
            a_ = 0;
            b_ = static_cast<double>(y_sum_);
            return;
        }

        if (static_cast<long double>(count_) * xx_sum_ - x_sum_ * x_sum_ == 0) {
            // all values in a bucket have the same key.
            a_ = 0;
            b_ = static_cast<double>(y_sum_) / count_;
            return;
        }

        auto slope = static_cast<double>(
                (static_cast<long double>(count_) * xy_sum_ - x_sum_ * y_sum_) /
                (static_cast<long double>(count_) * xx_sum_ - x_sum_ * x_sum_));
        auto intercept = static_cast<double>(
                (y_sum_ - static_cast<long double>(slope) * x_sum_) / count_);
        a_ = slope;
        b_ = intercept;

        // If floating point precision errors, fit spline
        if (a_ <= 0) {
            a_ = (y_max_ - y_min_) / (x_max_ - x_min_);
            b_ = -static_cast<double>(x_min_) * a_;
        }
	}

    void train(T *recs_, uint64_t size){
        T first_key = recs_[0];
        for (size_t i = size / 8; i < size * 7 / 8; i++){
            add(recs_[i] - first_key, i);
        }
        build();
    }

    double predict(double x) {
        return x * a_ + b_;
    }

private:
    int count_ = 0;
    long double x_sum_ = 0;
    long double y_sum_ = 0;
    long double xx_sum_ = 0;
    long double xy_sum_ = 0;
    T x_min_ = std::numeric_limits<T>::min();
    T x_max_ = std::numeric_limits<T>::max();
    double y_min_ = std::numeric_limits<double>::max();
    double y_max_ = std::numeric_limits<double>::lowest();
};

template <class T>
class WeightedLinearModelBuilder {
public:
    double a_;
    double b_;
    WeightedLinearModelBuilder() { }
    ~WeightedLinearModelBuilder() {delete[] weight;}
    inline void add(T x, int y, int w) {
        count_++;
        w_sum += w;
        x_sum_ += static_cast<long double>(w * x);
        y_sum_ += static_cast<long double>(w * y);
        xx_sum_ += static_cast<long double>(w * x) * x;
        xy_sum_ += static_cast<long double>(w * x) * y;
        x_min_ = std::min(x, x_min_);
        x_max_ = std::max(x, x_max_);
        y_min_ = std::min<double>(y, y_min_);
        y_max_ = std::max<double>(y, y_max_);
    }

    void build() {
        
        double slope = static_cast<double>(
                (static_cast<long double>(w_sum) * xy_sum_ - x_sum_ * y_sum_) /
                (static_cast<long double>(w_sum) * xx_sum_ - x_sum_ * x_sum_));
        double intercept = static_cast<double>(
                (y_sum_ - static_cast<long double>(slope) * x_sum_) / w_sum);
        a_ = slope;
        b_ = intercept;

        // If floating point precision errors, fit spline
        if (a_ <= 0) {
            a_ = (y_max_ - y_min_) / (x_max_ - x_min_);
            b_ = -static_cast<double>(x_min_) * a_;
        }
	}

    void cal_weight(T *recs_, uint64_t size){
        T first_key = recs_[0];
        int mid1_pos = (size - 1) / 3;
        int mid2_pos = (size - 1) * 2 / 3;

        const long double mid1_key =
            (static_cast<long double>(recs_[mid1_pos] - first_key) +
            static_cast<long double>(recs_[mid1_pos + 1] - first_key)) /
            2;
        const long double mid2_key =
            (static_cast<long double>(recs_[mid2_pos] - first_key) +
            static_cast<long double>(recs_[mid2_pos + 1] - first_key)) /
            2;
                
        const double mid1_target = mid1_pos;
        const double mid2_target = mid2_pos;

        long double slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
        long double intercept = mid1_target - slope * (mid1_key);

        // LinearModelBuilder<T> model;
        // model.train(recs_, size);
        // double slope = model.a_;
        // double intercept = model.b_;

        weight = new int[size];

        int i = 0, last_i = i;
        int cur_pos = predict(recs_[i], recs_[0], slope, intercept, size);
        for (i = 1; i < size; i++){
            int pre = predict(recs_[i], recs_[0], slope, intercept, size);
            if(pre != cur_pos){
                int count = i - last_i;
                for(int j = last_i; j < i; j++) weight[j] = count;
                last_i = i;
                cur_pos = pre;
            }
        }
        if(last_i != size){//the rest part
            int count = size - last_i;
            for(int j = last_i; j < i; j++) weight[j] = count;
        }
    }

    void train(T *recs_, uint64_t size){
        T first_key = recs_[0];
        cal_weight(recs_, size);
        for (size_t i = size / 8; i < size * 7 / 8; i++){
            add(recs_[i] - first_key, i, weight[i]);
        }
        build();
    }

    int predict(double x, const T &first_key, const double &slope, const double &intercept, const uint64_t &capacity) {
        double predict = (x - first_key) * slope + intercept + 0.5;
        if(likely(predict >= 0 && predict < capacity)) return predict;
        else if(predict >= capacity)    return capacity - 1;
        else return 0;
        // return x * a_ + b_;
    }

private:
    int count_ = 0;
    long double w_sum = 0;
    long double x_sum_ = 0;
    long double y_sum_ = 0;
    long double xx_sum_ = 0;
    long double xy_sum_ = 0;
    int *weight = nullptr;
    T x_min_ = std::numeric_limits<T>::min();
    T x_max_ = std::numeric_limits<T>::max();
    double y_min_ = std::numeric_limits<double>::max();
    double y_max_ = std::numeric_limits<double>::lowest();
};

template<class T>
class FMCD
{
public:
    double a_;
    double b_;
    int D = 1;
    bool success = false;
public:
    FMCD(){
        a_ = 0;
        b_ = 0;
    }
    void train(T *recs_, uint64_t size){
        // FMCD method
        // Here the implementation is a little different with Algorithm 1 in our
        // paper. In Algorithm 1, U_T should be (keys[size-1-D] - keys[D]) / (L
        // - 2). But according to the derivation described in our paper, M.A
        // should be less than 1 / U_T. So we added a small number (1e-6) to
        // U_T. In fact, it has only a negligible impact of the performance.
        T first_key = recs_[0];
        success = false;
            {
            const int L = size;
            int i = 0;
            D = 1;
            double Ut = (static_cast<long double>(recs_[size - 1 - D] - first_key) -
                        static_cast<long double>(recs_[D] - first_key)) /
                            (static_cast<double>(L - 2)) +
                        1e-6;
            while (i < size - 1 - D) {
                while (i + D < size && recs_[i + D] - recs_[i] >= Ut) {
                    i++;
                }
                if (i + D >= size) {
                    break;
                }
                D = D + 1;
                if (D * 3 > size)
                break;
                Ut = (static_cast<long double>(recs_[size - 1 - D] - first_key) -
                    static_cast<long double>(recs_[D] - first_key)) /
                        (static_cast<double>(L - 2)) +
                    1e-6;
            }
            if (D * 3 <= size) {

                a_ = 1.0 / Ut;
                b_ =
                    (L -
                    a_ * (static_cast<long double>(recs_[size - 1 - D] - first_key) +
                                    static_cast<long double>(recs_[D] - first_key))) /
                    2;
                fmcd_success++;
                fmcd_D += D;
                success = true;
            } else {

                int mid1_pos = (size - 1) / 3;
                int mid2_pos = (size - 1) * 2 / 3;


                const long double mid1_key =
                    (static_cast<long double>(recs_[mid1_pos] - first_key) +
                    static_cast<long double>(recs_[mid1_pos + 1] - first_key)) /
                    2;
                const long double mid2_key =
                    (static_cast<long double>(recs_[mid2_pos] - first_key) +
                    static_cast<long double>(recs_[mid2_pos + 1] - first_key)) /
                    2;

                
                const double mid1_target =
                    mid1_pos;
                const double mid2_target =
                    mid2_pos;

                a_ = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                b_ = mid1_target - a_ * (mid1_key);
                fmcd_failure++;
            }
            }
    }
    ~FMCD(){}
};

struct segment{
    uint64_t start = 0;
    uint64_t end = 0;
    segment(uint64_t s, uint64_t e){
        start = s;
        end = e;
    }
    inline uint64_t number(){
        return end - start + 1;
    }    
};
template<class T, class P>
class SecondDerivative{//get segment using second derivative
public:

public:
    std::vector<segment> segments;
    SecondDerivative(T *recs_, uint64_t size, uint64_t num_seg = 20000)  {
        if(num_seg == 1){
            segments.emplace_back(segment(0, size - 1));
            return;
        }
        segments.reserve(num_seg);
        auto cmp = [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
            return a.first > b.first;
        };
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(cmp)> heap(cmp);
        T second_derivative;
        int no_seg_range = 8;//size / 4;
        for (size_t i = no_seg_range; i < size - no_seg_range; i++){//calculate the second derivation of each pos
            second_derivative = labs((T)(recs_[i] - recs_[i-1]) - (recs_[i-1] - recs_[i-2]));
            if(heap.size() == num_seg - 1 && second_derivative > heap.top().first){
                heap.pop();
                heap.emplace(std::make_pair(second_derivative, i - 1));
            }
            else if(heap.size() < num_seg -1 ){
                heap.emplace(std::make_pair(second_derivative, i - 1));
            }
        }

        std::vector<uint8_t> split_index(size, 0);
        split_index[size - 1] = 1;
        
        while(!heap.empty()){
            split_index[heap.top().second] = 1;
            heap.pop();
        }

        int start_i = 0;
        for (size_t i = 0; i < size; i++){
            if(split_index[i] == 1){
                // while( i < size - 1 && recs_[i] == recs_[i + 1]) i++;//handle the duplication
                if(i - start_i + 1 < 8){//the segment should not be too small
                    if(start_i == 0)
                        continue;
                    segments.back().end = i;
                }
                else{
                    segments.emplace_back(segment(start_i, i));
                }
                start_i = i + 1;
            }
        }
    }
    ~SecondDerivative(){}

    inline uint64_t number(){
        return segments.size();
    }

    segment operator [](int i){
        return segments[i];
    }
};
//探测该数据段是否会导致溢出数量查过一定阈值，如果超过则需要对数据进行进一步分割
#define SPLIT_THRESHOLD (12<<6)
#define MERGE_THRESHOLD (1<<4)
static void MergeSeg(std::vector<uint64_t> &segments){//合并元素数量过少的数据段
    size_t size = segments.size();
    if(size <= 1) return;
    for (size_t i = 0; i < size; i++){
        if(i >= 0 && i < size - 1){
            if(segments[i] < MERGE_THRESHOLD){
                segments[i+1] += segments[i];
                segments[i] = 0;
            }
        }
        else if(i == size - 1){
            if(segments[i] < MERGE_THRESHOLD){
                segments[i-1] += segments[i];
                segments[i] = 0;
            }
        }
    }
    
}
struct NodeParamer
{
    uint32_t size;
    uint32_t capacity;
    double a_, b_;
    NodeParamer(uint32_t s, uint32_t c, double a, double b):
    size(s), capacity(c), a_(a), b_(b){};
    NodeParamer(){}
};

template<class T>
__always_inline int Predict(const T &key, const T &first_key, const double &slope, const double &intercept, const uint64_t &capacity){
    double predict = (key - first_key) * slope + intercept + 0.5;
    if(likely(predict >= 0 && predict < capacity)) return predict;
    else if(predict >= capacity)    return capacity - 1;
    else return 0;
    __builtin_unreachable();

}
template<class T>
using LeafModel = WeightedLinearModelBuilder<T>;
template<class T>
using InnerModel = LinearModelBuilder<T>;

template<class T, class P>
std::vector<NodeParamer> TrySplit(T *recs_, uint64_t size){
    std::vector<NodeParamer> result;
    result.reserve(10);
    LeafModel<T> model;
    // FMCD<T> model;
    model.train(recs_, size);
    uint64_t capacity = size * SCALE_FACTOR ;
    double slope = model.a_ * SCALE_FACTOR;
    double intercept = model.b_ *SCALE_FACTOR;

    int i = 0, last_i = i;
    int cur_pos = Predict(recs_[i], recs_[0], slope, intercept, capacity);
    for (i = 1; i < size; i++){
        int predict = Predict(recs_[i], recs_[0], slope, intercept, capacity);
        if(predict != cur_pos){
            int count = i - last_i;
            if(count >= SPLIT_THRESHOLD){
                SecondDerivative<T, P> segments(recs_, size, 2);
                uint64_t start_pos = 0;
                for (size_t j = 0; j < segments.number(); j++){
                    std::vector<NodeParamer> sub_result = TrySplit<T, P>(&recs_[start_pos], segments[j].number());
                    result.insert(result.end(), sub_result.begin(), sub_result.end());
                    start_pos += segments[j].number();
                }
                return result;
            }
            last_i = i;
            cur_pos = predict;
        }
    }

    if(last_i != size){//the rest part
        int count = size - last_i;
        if(count >= SPLIT_THRESHOLD){
            SecondDerivative<T, P> segments(recs_, size, 2);
            uint64_t start_pos = 0;
            for (size_t i = 0; i < segments.number(); i++){
                std::vector<NodeParamer> sub_result = TrySplit<T, P>(&recs_[start_pos], segments[i].number());
                result.insert(result.end(), sub_result.begin(), sub_result.end());
                start_pos += segments[i].number();
            }
            return result;
        }
    }
    result.emplace_back(NodeParamer(size, capacity, slope, intercept));
    return result;
}

template<class T, class P>
int BinarySearch(Record<T,P> *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        #ifndef NDEBUG
        of_probe_length++;
        #endif
        mid = low + (high - low) / 2;
        if(recs_[mid].key < target_key) low = mid + 1;
        else if (recs_[mid].key > target_key) high = mid;
        else return mid;
    }
    return -1;
    
};

//find the first position greater than target_key
template<class T, class P>
int BinarySearch_UpperBound(Record<T,P> *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        mid = low + (high - low) / 2;
        if(recs_[mid].key <= target_key) low = mid + 1;
        else if (recs_[mid].key > target_key) high = mid;
        // else return mid;
    }
    return low;    
};
//find the first position greater than or equal to target_key
template<class T, class P>
int BinarySearch_LowerBound(Record<T,P> *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        mid = low + (high - low) / 2;
        if(recs_[mid].key < target_key) low = mid + 1;
        else if (recs_[mid].key >= target_key) high = mid;
        // else return mid;
    }
    return low;
};
//branchless binary search
template <class T, class P>
int BinarySearch_Branchless(Record<T,P> *recs_, int size, T target_key)
{
    int left = 0, right = (int)size - 1; // both inclusive
    while (left <= right)
    {
        int mid = left + (right - left) / 2;
        if (target_key == recs_[mid].key) return mid;
        const int midMinusOne = mid - 1;
        const int midPlusOne = mid + 1;
        right = target_key <= recs_[mid].key ? midMinusOne : right;
        left = target_key <= recs_[mid].key ? left : midPlusOne;
    }
    return -1; // not found
}

template <class T, class P>
int BinarySearch_LowerBound_Branchless(Record<T,P> *recs_, int size, T target_key)
{
    int left = 0, right = size; // left inclusive, right exclusive
    while (left < right)
    {

        int mid = left + (right - left) / 2;
        const int midPlusOne = mid + 1;
        right = target_key <= recs_[mid].key ? mid : right;
        left = target_key <= recs_[mid].key ? left : midPlusOne;
    }
    return left;
}

template <class T, class P>
int BinarySearch_UpperBound_Branchless(Record<T,P> *recs_, int size, T target_key)
{
    int left = 0, right = size; // left inclusive, right exclusive
    while (left < right)
    {

        int mid = left + (right - left) / 2;
        const int midPlusOne = mid + 1;
        right = target_key < recs_[mid].key ? mid : right;
        left = target_key < recs_[mid].key ? left : midPlusOne;
    }
    return left;
}

template<class T>
int BinarySearch(T *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        #ifndef NDEBUG
        of_probe_length++;
        #endif
        mid = low + (high - low) / 2;
        if(recs_[mid] < target_key) low = mid + 1;
        else if (recs_[mid] > target_key) high = mid;
        else return mid;
    }
    return -1;
    
};

//find the first position greater than target_key
template<class T>
int BinarySearch_UpperBound(T *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        mid = low + (high - low) / 2;
        if(recs_[mid] <= target_key) low = mid + 1;
        else if (recs_[mid] > target_key) high = mid;
        // else return mid;
    }
    return low;    
};
//find the first position greater than or equal to target_key
template<class T>
int BinarySearch_LowerBound(T *recs_, int length, T target_key){
    int low = 0;
    int high = length;
    int mid;
    while (low < high){
        mid = low + (high - low) / 2;
        if(recs_[mid] < target_key) low = mid + 1;
        else if (recs_[mid] >= target_key) high = mid;
        // else return mid;
    }
    return low;
};


}//namespace lisk
#endif