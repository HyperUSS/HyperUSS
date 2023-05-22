# CPU Code

## Repository Structure

- `Common/`: hash and mmap functions and basic settings for experiments
- `Struct/`: data structures for USS
- `Multiple/`: the implemented version of HyperUSS, Single Cocosketch and USS
- `BenchMark.h`: benchmark for experiments on real dataset
- `BenchMark_HeavyChange.h`: benchmark for heavychange experiments
- `BenchMark_Synthetic.h`: benchmark for experiments on synthetic datasets
- `PartialKey.txt`: a text file for subset queries
- `main.cpp`: an example on a 50m_10K10V synthetic dataset using 500KB memory

## How to run

First of all, generate a synthetic dataset (50m items, 10 integers for key and 10 attributes) in `data/Synthetic/`.

```
mkdir build
cd build
cmake ..
make
./CPU ../../data/Synthetic/Hyper_zipf_1.500_50m_10K10V.txt
```

Results can be found in `example_FullKey.txt` and `example_PartialKey.txt`.
