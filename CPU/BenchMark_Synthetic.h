#ifndef HHBENCHSYN_H
#define HHBENCHSYN_H

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "Hyper.h"
#include "MMap.h"
#include "Single.h"
#include "USS.h"

class BenchMark_Synthetic {
 public:
  BenchMark_Synthetic(std::string PATH, std::string FILENAME,
                      double SKEW_FACTOR = 1) {
    filename = FILENAME;
    skew_factor = SKEW_FACTOR;
    std::ifstream input(PATH);
    if (!input) {
      std::cout << "fail to open file " << PATH << std::endl;
    }
    while (!input.eof()) {
      TUPLES tuple;
#ifndef AVERAGE
      for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        uint32_t tmp;
        input >> std::dec >> tmp;
        tuple.value.values[i] = tmp;
#ifdef SKEW
        if (i == 0) {
          tuple.value.values[i] *= SKEW_FACTOR;
        }
#endif
      }
#else
      tuple.value.values[0] = 1.0;
      for (int i = 1; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        uint32_t tmp;
        input >> std::dec >> tmp;
        tuple.value.values[i] = tmp;
      }
#endif
      if (input.eof()) {
        break;
      }
      for (int j = 0; j < TUPLES_NUM; ++j) {
        input >> std::hex >> tuple.id.key[j];
      }
      dataset.emplace_back(tuple);
    }

    std::cout << "read complete!" << std::endl;
    length = dataset.size();

    for (int i = 0; i < length; ++i) {
      tuplesMp[dataset[i].id] += dataset[i].value;
    }

    std::cout << "Number of different keys: " << tuplesMp.size() << std::endl;

    uint32_t key_num = tuplesMp.size();
    for (auto &[key, value] : tuplesMp) {
      key_vector.push_back(key);
    }

    srand(time(NULL));
    for (int i = 0; i < SUBSETS_NUM; ++i) {
      std::set<uint32_t> tmp;
      //   while (tmp.size() < SUBSETS_SIZE_ALPHA * key_num) {
      while (tmp.size() < 1) {
        uint32_t rd = rand() % key_num;
        tmp.insert(rd);
      }
      subsets.push_back(tmp);
    }
  }

  ~BenchMark_Synthetic() {}

  void HHHyperBench(uint32_t MEMORY, uint32_t HASH_NUM = 2,
                    std::string BENCHNAME = "Ours") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start HyperBench!" << std::endl;

    HASH_NUM = 2;

    OurHyper *sketch = new OurHyper(MEMORY, HASH_NUM);

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (uint32_t i = 0; i < length; ++i) {
      sketch->Insert(dataset[i]);
    }

    auto end_insert = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_insert =
        end_insert - start_insert;
    throughput = length / (tm_insert.count() / 1000);
    std::cout << "Hyper Insert Time: " << tm_insert.count() / 1000 << "s"
              << std::endl;
    std::cout << "Hyper Throughput: " << throughput << std::endl;

    // query
    auto start_query = std::chrono::high_resolution_clock::now();
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple = sketch->AllQuery();
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "Hyper Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

#ifndef SKEW
    CompareSubset(estTuple, tuplesMp);
#else
    CompareSubset_SKEW(estTuple, tuplesMp);
#endif

    delete sketch;
  }

  void HHSingleBench(uint32_t MEMORY, uint32_t HASH_NUM = 2,
                     std::string BENCHNAME = "Coco") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start SingleBench!" << std::endl;

    OurSingle *sketch[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch[i] = new OurSingle(MEMORY / TUPLES_VALUES_ELEMENT_NUM, HASH_NUM);
    }

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < length; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        Single tmp;
        tmp.id = dataset[i].id;
        tmp.value = dataset[i].value.values[j];
        sketch[j]->Insert(tmp);
      }
    }
    auto end_insert = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_insert =
        end_insert - start_insert;
    throughput = length / (tm_insert.count() / 1000);
    std::cout << "Single Insert Time: " << tm_insert.count() / 1000 << "s"
              << std::endl;
    std::cout << "Single Throughput: " << throughput << std::endl;

    // query
    auto start_query = std::chrono::high_resolution_clock::now();
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple;
    std::vector<std::unordered_map<TUPLES_ID, double> > tmpTuple(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple[i] = sketch[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple[i].begin(); it != tmpTuple[i].end(); ++it) {
        estTuple[it->first].values[i] = it->second;
      }
    }
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "Single Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

    CompareSubset(estTuple, tuplesMp);

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      delete sketch[i];
    }
  }

  void HHUSSBench(uint32_t MEMORY, std::string BENCHNAME = "USS") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start USSBench!" << std::endl;

    OurUSS *sketch[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch[i] = new OurUSS(MEMORY / TUPLES_VALUES_ELEMENT_NUM);
    }

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < length; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        sketch[j]->Insert(dataset[i].id, dataset[i].value.values[j]);
      }
    }
    auto end_insert = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_insert =
        end_insert - start_insert;
    throughput = length / (tm_insert.count() / 1000);
    std::cout << "USS Insert Time: " << tm_insert.count() / 1000 << "s"
              << std::endl;
    std::cout << "USS Throughput: " << throughput << std::endl;

    // query
    auto start_query = std::chrono::high_resolution_clock::now();
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple;
    std::vector<std::unordered_map<TUPLES_ID, double> > tmpTuple(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple[i] = sketch[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple[i].begin(); it != tmpTuple[i].end(); ++it) {
        estTuple[it->first].values[i] = it->second;
      }
    }
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "USS Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

    CompareSubset(estTuple, tuplesMp);

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      delete sketch[i];
    }
  }

 private:
  std::string filename;
  std::string benchname;

  std::vector<TUPLES> dataset;

  double skew_factor;
  uint32_t length;
  uint32_t memory;
  double throughput;

  std::vector<TUPLES_ID> key_vector;
  std::vector<std::set<uint32_t> > subsets;

  std::unordered_map<TUPLES_ID, TUPLES_VALUE> tuplesMp;

  template <class T>
  // mp : estimate
  // record : real
  void CompareSubset(T mp, T record) {
    double aae_sum = 0.0;
    double are_sum = 0.0;
#ifndef AVERAGE
    double NUM = 0;
    for (int i = 0; i < SUBSETS_NUM; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        double mp_sum = 0.0, record_sum = 0.0;
        for (auto x : subsets[i]) {
          mp_sum += mp[key_vector[x]].values[j];
          record_sum += record[key_vector[x]].values[j];
        }
        if (record_sum > 0) {
          NUM += 1;
          aae_sum += abs(mp_sum - record_sum);
          are_sum += abs(mp_sum - record_sum) / record_sum;
        }
      }
    }
#else
    double NUM = 0;
    std::vector<double> mp_count(SUBSETS_NUM, 0), record_count(SUBSETS_NUM, 0);
    for (int i = 0; i < SUBSETS_NUM; ++i) {
      for (auto x : subsets[i]) {
        mp_count[i] += mp[key_vector[x]].values[0];
        record_count[i] += record[key_vector[x]].values[0];
      }
    }
    for (int i = 0; i < SUBSETS_NUM; ++i) {
      for (int j = 1; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        double mp_sum = 0.0, record_sum = 0.0;
        for (auto x : subsets[i]) {
          mp_sum += mp[key_vector[x]].values[j];
          record_sum += record[key_vector[x]].values[j];
        }
        if (mp_count[i] > 0 && record_count[i] > 0 && record_sum > 0) {
          NUM += 1;
          double mp_avg = mp_sum / mp_count[i];
          double record_avg = record_sum / record_count[i];
          aae_sum += abs(mp_avg - record_avg);
          are_sum += abs(mp_avg - record_avg) / record_avg;
        }
      }
    }
#endif
    std::ofstream ofs;
    ofs.open(filename + ".txt", std::ios::app);
    ofs << benchname << " " << memory << " " << aae_sum / NUM << " "
        << are_sum / NUM << " " << throughput << std::endl;
    ofs.close();
  }

  // SKEW START HERE
  template <class T>
  // mp : estimate
  // record : real
  void CompareSubset_SKEW(T mp, T record) {
    std::vector<double> are_vec(TUPLES_VALUES_ELEMENT_NUM, 0);
    std::vector<double> cnt_vec(TUPLES_VALUES_ELEMENT_NUM, 0);
    for (int i = 0; i < SUBSETS_NUM; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        double mp_sum = 0.0, record_sum = 0.0;
        for (auto x : subsets[i]) {
          mp_sum += mp[key_vector[x]].values[j];
          record_sum += record[key_vector[x]].values[j];
        }
        if (record_sum > 0) {
          cnt_vec[j] += 1;
          are_vec[j] += abs(mp_sum - record_sum) / record_sum;
        }
      }
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      are_vec[i] /= cnt_vec[i];
    }
    std::ofstream ofs;
    ofs.open(filename + ".txt", std::ios::app);
    ofs << memory << " " << skew_factor << " "
        << *std::min_element(are_vec.begin(), are_vec.end()) << " "
        << *std::max_element(are_vec.begin(), are_vec.end()) << " "
        << std::accumulate(are_vec.begin(), are_vec.end(), 0.0) /
               TUPLES_VALUES_ELEMENT_NUM
        << std::endl;
    ofs.close();
  }
};

#endif
