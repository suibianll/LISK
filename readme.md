
Implementation of the paper "LISK: A High-Performance Learned Index for String Keys".

### Compling

Assuming to compile under a `build` directory:
```bash
git clone https://github.com/baotonglu/dex.git
cd LISK
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. 
make -j
```
### Running test

```bash
cd build
#test lookup
./microbench dataset
#test insert
./loadbench dataset
```