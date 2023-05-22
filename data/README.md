# Zipf Dataset

We generate a series of synthetic datasets for experiment. The keys in the synthetic dataset are randomly generated integers, following the Zipf distribution with skewness of 1.5. For each item, the value of each attribute is a randomly positive integer, following a certain exponential distribution.

## How to run

```
cd Synthetic
g++ main.cpp -o main && ./main
```

You can change parameters of the synthetic dataset by modifying `main.cpp` correspondingly.
