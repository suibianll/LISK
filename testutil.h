
#ifndef __TEST_UTIL__
#define __TEST_UTIL__

#include <sstream>
#include <iostream>
#include <functional>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <signal.h>
#include "zipf.h"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <vector>

template<class T>
long long load_binary_data(T *&data, long long length, const std::string &file_path) {
    // open key file
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if (!is.is_open()) {
        return 0;
    }

    std::cout << file_path << std::endl;

    // read the number of keys
    T max_size;
    is.read(reinterpret_cast<char*>(&max_size), sizeof(T));

    std::cout << max_size << std::endl;

    // create array
    if(length < 0 || length > max_size) length = max_size;
    data = new T[length];

    // read keys
    is.read(reinterpret_cast<char *>(data), std::streamsize(length * sizeof(T)));
    is.close();
    return length;
}

template<class T>
long long load_text_data(T *&array, long long length, const std::string &file_path) {
    std::ifstream is(file_path.c_str());
    if (!is.is_open()) {
        return 0;
    }
    long long i = 0;
    std::string str;

    std::vector<T> temp_keys;
    temp_keys.reserve(200000000);
    while (std::getline(is, str) && (i < length || length < 0)) {
        std::istringstream ss(str);
        T key;
        ss >> key;
        temp_keys.push_back(key);
        i++;
    }

    array = new T[temp_keys.size()];
    for(int j = 0; j < temp_keys.size(); j++) {
        array[j] = temp_keys[j];
    }
    is.close();
    return temp_keys.size();
}

template<class T>
T *get_search_keys(T array[], int num_keys, int num_searches, size_t *seed = nullptr) {
    auto *keys = new T[num_searches];

#pragma omp parallel
    {
        std::mt19937_64 gen(std::random_device{}());
        if (seed) {
            gen.seed(*seed + omp_get_thread_num());
        }
        std::uniform_int_distribution<int> dis(0, num_keys - 1);
#pragma omp for
        for (int i = 0; i < num_searches; i++) {
            int pos = dis(gen);
            keys[i] = array[pos];
        }
    }

    return keys;
}


bool file_exists(const std::string &str) {
    std::ifstream fs(str);
    return fs.is_open();
}

template<class T>
T *get_search_keys_zipf(T array[], int num_keys, int num_searches, size_t *seed = nullptr) {
    auto *keys = new T[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys, seed);
    for (int i = 0; i < num_searches; i++) {
        int pos = zipf_gen.nextValue();
        keys[i] = array[pos];
    }
    return keys;
}


template<typename T>
T *unique_data(T *key1, size_t &size1, T *key2, size_t &size2) {
    size_t ptr1 = 0;
    size_t ptr2 = 0;

    std::sort(key1, key1 + size1);
    size1 = std::unique(key1, key1 + size1) - key1;
    std::sort(key2, key2 + size2);
    size2 = std::unique(key2, key2 + size2) - key2;

    size_t result = 0;

    while (ptr1 < size1 && ptr2 < size2) {
        while (key1[ptr1] < key2[ptr2] && ptr1 < size1) {
            ptr1++;
        }
        if (key1[ptr1] == key2[ptr2]) {
            ptr2++;
            continue;
        }
        key2[result++] = key2[ptr2++];
    }

    while (ptr2 < size2) {
        key2[result++] = key2[ptr2++];
    }

    size2 = result;
    std::random_shuffle(key2, key2 + size2);

    return &key2[result];
}
inline double get_now() { 
    struct timeval tv; 
    gettimeofday(&tv, 0); 
    return tv.tv_sec + tv.tv_usec / 1000000.0; 
} 

#endif