#include <bits/stdc++.h>

#include "generator.h"
#include "murmur3.h"

using namespace std;

uint32_t key_num = 5;
uint32_t value_num = 2;

uint32_t key_len = 4;
uint32_t value_len = 4;

uint32_t flow_num = 5000;
uint64_t packet_num = 5E7;

vector<double> value_lambda = {1, 2, 4, 8, 16, 1, 2, 4, 8, 16};
// vector<double> value_lambda = {128, 128};

void gen_zipf_dataset(double alpha) {
  string filename =
      "Hyper_zipf_" + to_string(alpha).substr(0, 5) + "_10K10V.txt";
  ofstream outFile(filename.c_str());

  std::cout << filename << " start" << std::endl;

  for (uint64_t i = 0; i < packet_num; ++i) {
    // value: random value
    for (int j = 0; j < value_num; ++j) {
      uint32_t rand_value = (uint32_t)(exponential(1.0 / value_lambda[j % 5]));
      outFile << std::dec << rand_value << '\t';
    }

    // key
    uint32_t rand_num = zipf(alpha, flow_num, i == 0);
    for (int j = 1; j <= key_num; ++j) {
      uint32_t key = MurmurHash3_x86_32(&rand_num, key_len, j);
      outFile << std::hex << key << '\t';
    }

    outFile << '\n';
  }

  std::cout << filename << " done" << std::endl;
}

void gen_skew_zipf_dataset(double alpha) {
  string filename =
      "Hyper_skew_zipf_" + to_string(alpha).substr(0, 5) + "_5K5V.txt";
  ofstream outFile(filename.c_str());

  std::cout << filename << " begin" << std::endl;

  for (uint64_t i = 0; i < packet_num; ++i) {
    // value: 1
    uint32_t value = 1;
    outFile << std::dec << value << '\t';

    // value: [0, 100)
    uint32_t rand_value1 = rand() % 100;
    outFile << std::dec << rand_value1 << '\t';

    // value: [0, 10000)
    uint32_t rand_value2 = rand() % 10000;
    outFile << std::dec << rand_value2 << '\t';

    // value: [0, 1000000000)
    uint32_t rand_value3 = rand() % 1000000000;
    outFile << std::dec << rand_value3 << '\t';

    // value: [1000000000, 2000000000)
    uint32_t rand_value4 = rand() % 1000000000 + 1000000000;
    outFile << std::dec << rand_value4 << '\t';

    // key
    uint32_t rand_num = zipf(alpha, flow_num, i == 0);
    for (int j = 1; j <= key_num; ++j) {
      uint32_t key = MurmurHash3_x86_32(&rand_num, key_len, j);
      outFile << std::hex << key << '\t';
    }

    outFile << '\n';
  }

  std::cout << filename << " done" << std::endl;
}

int main() {
  //   gen_zipf_dataset(1.0);

  gen_zipf_dataset(1.5);

  //   gen_skew_zipf_dataset(1.0);

  //   gen_skew_zipf_dataset(1.5);
}
