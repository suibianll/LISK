
#include <iostream>
#include <fstream>
#include "testutil.h"
#include "./src/lisk.h"


using namespace lisk;
size_t table_size;
size_t load_size;

int main(int argc, char* argv[]){
    LISK *index = new LISK();
    std::vector<Key> keys;
    std::vector<_value_t> values;
    std::string file = argv[1];
    std::string filepath = "/path/to/" + file;
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

    uint64_t opnum = 200000000;
    Key* query_keys = get_search_keys_zipf(keys.data(), load_size, opnum);
    // _key_t* query_keys = get_search_keys(init_data, table_size, opnum);

    double start_time;
    double end_time;
    uint64_t val;
    uint64_t failed = 0;
    start_time = get_now();

    for (size_t i = 0; i < opnum; i++)
    {
        bool found = index->find(query_keys[i], val);
        #ifndef NDEBUG
        if(!found) {
            ++failed;
            found = index->find(query_keys[i], val);
        }
        #endif
    }

    end_time = get_now();
    double tput = opnum / (end_time - start_time) / 1000000; //Mops/sec
    std::cout << file << ": throughput " << tput << " Mops/s"<<std::endl;
    
    #ifndef NDEBUG
    std::cout<<"average search depth: "<<(double)lisk::index_depth/opnum<<std::endl;
    std::cout<<"average inner probe length: "<<(double)lisk::inner_probe_length/opnum<<std::endl;
    std::cout<<"average leaf probe length: "<<(double)lisk::leaf_probe_length/opnum<<std::endl;
    std::cout<<"average overflow probe length: "<<(double)lisk::of_probe_length/lisk::of_search_time<<std::endl;
    std::cout<<"total overflow search time: "<<lisk::of_search_time<<std::endl;
    std::cout<<"failed search: "<<failed<<std::endl;
    #endif
    #ifndef NSTAT
    std::cout << "latency breakdown(us)"<<std::endl;
    std::cout << "inner search:"<< lisk::latency_breakdown[lisk::INNERTIME]/opnum<<std::endl;
    std::cout << "leaf search:"<< lisk::latency_breakdown[lisk::LEAFTIME]/opnum<<std::endl;
    std::cout << "leaf probe:"<< lisk::latency_breakdown[lisk::LEAFPROBE]/opnum<<std::endl;
    std::cout << "overflow search:"<< lisk::latency_breakdown[lisk::OVERFLOWPROBE]/lisk::of_search_time<<std::endl;
    #endif

}