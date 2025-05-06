
#include <iostream>
#include <fstream>
#include <unordered_set>
#include "testutil.h"
#include "./src/lisk.h"

using namespace lisk;

void loadbench(std::string file){
    google::InitGoogleLogging("loadtest");
    google::SetLogDestination(google::GLOG_INFO, "/home/czl/trie_index/sbench/LISK/log/INFO-");
    google::SetLogDestination(google::GLOG_ERROR, "/home/czl/trie_index/sbench/LISK/log/ERROR-");
    google::SetStderrLogging(google::GLOG_FATAL);
    std::string init_file = "/home/czl/trie_index/sbench/ycsb/workloads/" + file + "_load.dat";
    std::string txn_file = "/home/czl/trie_index/sbench/ycsb/workloads/" + file + "_query_l.dat";
    std::vector<Key> init_keys;
    std::vector<_value_t> init_values;
    std::vector<Key> query_keys;
    std::ifstream infile_load(init_file);
    std::string op;
    std::string key_str;
    std::string line;
    std::unordered_set<std::string> dedup;
    dedup.reserve(10000000);
    while (getline(infile_load, line)){
        if(line.empty()) continue;
        std::istringstream iss(line);
        Key k;
        iss >> op >> key_str;
        if(dedup.count(key_str) != 0) continue;
        else dedup.insert(key_str);
        k.set(key_str.c_str(), key_str.size());
        init_keys.push_back(k);
    }

    std::ifstream infile_txn(txn_file);
    while (getline(infile_txn, line)){
        if(line.empty()) continue;
        std::istringstream iss(line);
        Key k;
        iss >> op >> key_str;
        k.set(key_str.c_str(), key_str.size());
        query_keys.push_back(k);
    }
    std::cout<<"load size:"<<init_keys.size()<<" insert size:"<<query_keys.size()<<std::endl;
    LISK *index = new LISK();
    std::sort(init_keys.begin(), init_keys.end());
    for (size_t i = 0; i < init_keys.size(); i++){
        init_values.emplace_back((_value_t)(&init_keys[i]));
        init_keys[i].setVal(init_values.back());
    }
    index->bulk_load(init_keys, init_values, init_keys.size());

    uint64_t opnum = query_keys.size();

    double start_time;
    double end_time;
    uint64_t val;
    start_time = get_now();

    for (size_t i = 0; i < opnum; i++){
        _value_t v = (_value_t)(&query_keys[i]);
        index->insert(query_keys[i], v);
    }

    end_time = get_now();
    double tput = opnum / (end_time - start_time) / 1000000; //Mops/sec
    FILE *fp=fopen("/home/czl/trie_index/sbench/LISK/test/statistics_load.txt","a+");
    std::cout << file << ": throughput " << tput << " Mops/s"<<std::endl;
    fprintf(fp, "%s: throughput %lf Mops/s\n",file.c_str(), tput);
}

int main(int argc, char* argv[]){
    std::string file = argv[1];
    // loadbench(file);
    LISK *index = new LISK();
    std::vector<Key> keys;
    std::vector<_value_t> values;
    std::string filepath = "/home/czl/dataset/vary_len/" + file;
    keys.reserve(10000000);
    values.reserve(10000000);
    std::ifstream in(filepath, std::ios::in);
    std::string line;
    while (std::getline(in, line)){
        if(line.empty()) continue;
        Key k;
        k.set(line.c_str(), line.size());
        keys.emplace_back(k);
    }
    in.close();
        
    int load_size = keys.size() * 0.5;
    std::sort(keys.begin(), keys.begin() + load_size);
    for (size_t i = 0; i < keys.size(); i++){
        values.emplace_back((_value_t)(&keys[i]));
        keys[i].setVal(values.back());
    }

    index->bulk_load(keys, values, load_size);

    uint64_t opnum = keys.size() - load_size;

    double start_time;
    double end_time;
    uint64_t val;
    uint64_t failed = 0;
    start_time = get_now();

    for (size_t i = 0; i < opnum; i++){
        index->insert(keys[load_size + i], values[load_size + i]);
    }

    end_time = get_now();
    double tput = opnum / (end_time - start_time) / 1000000; //Mops/sec
    FILE *fp=fopen("/home/czl/trie_index/sbench/LISK/test/statistics_load.txt","a+");
    std::cout << file << ": throughput " << tput << " Mops/s"<<std::endl;
    fprintf(fp, "%s: throughput %lf Mops/s\n",file.c_str(), tput);
    #ifndef NDEBUG
    std::cout<<"average search depth: "<<(double)lisk::index_depth/opnum<<std::endl;
    std::cout<<"average inner probe length: "<<(double)lisk::inner_probe_length/opnum<<std::endl;
    std::cout<<"average leaf probe length: "<<(double)lisk::leaf_probe_length/opnum<<std::endl;
    std::cout<<"average overflow probe length: "<<(double)lisk::of_probe_length/lisk::of_search_time<<std::endl;
    std::cout<<"total overflow search time: "<<lisk::of_search_time<<std::endl;
    std::cout<<"failed search: "<<failed<<std::endl;
    fprintf(fp, "average search depth:  %lf\n", (double)lisk::index_depth/opnum);
    fprintf(fp, "average inner probe length:  %lf\n", (double)lisk::inner_probe_length/opnum);
    fprintf(fp, "average leaf probe length:  %lf\n", (double)lisk::leaf_probe_length/opnum);
    fprintf(fp, "average overflow probe length:  %lf\n", (double)lisk::of_probe_length/lisk::of_search_time);
    fprintf(fp, "total overflow search time:  %ld\n", lisk::of_search_time);
    #endif
    #ifndef NSTAT
    std::cout << "latency breakdown(us)"<<std::endl;
    std::cout << "inner search:"<< lisk::latency_breakdown[lisk::INNERTIME]/opnum<<std::endl;
    std::cout << "leaf search:"<< lisk::latency_breakdown[lisk::LEAFTIME]/opnum<<std::endl;
    std::cout << "leaf probe:"<< lisk::latency_breakdown[lisk::LEAFPROBE]/opnum<<std::endl;
    std::cout << "overflow search:"<< lisk::latency_breakdown[lisk::OVERFLOWPROBE]/lisk::of_search_time<<std::endl;
    fprintf(fp, "latency breakdown(us)\n");
    fprintf(fp, "inner search:  %lf\n", lisk::latency_breakdown[lisk::INNERTIME]/opnum);
    fprintf(fp, "leaf search:  %lf\n", lisk::latency_breakdown[lisk::LEAFTIME]/opnum);
    fprintf(fp, "leaf probe:  %lf\n", lisk::latency_breakdown[lisk::LEAFPROBE]/opnum);
    fprintf(fp, "overflow search:  %lf\n", lisk::latency_breakdown[lisk::OVERFLOWPROBE]/lisk::of_search_time);
    #endif

}