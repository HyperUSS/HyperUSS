#ifndef HHBENCHHC_H
#define HHBENCHHC_H

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

class BenchMark_HeavyChange {
 public:
  BenchMark_HeavyChange(std::string PATH, std::string FILENAME,
                        double ALPHA = 0.000001) {
    alpha = ALPHA;
    filename = FILENAME;
    std::ifstream input(PATH);
    if (!input) {
      std::cout << "fail to open file " << PATH << std::endl;
    }
    while (!input.eof()) {
      TUPLES tuple;
      for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        double tmp;
        input >> std::dec >> tmp;
        tuple.value.values[i] = tmp;
      }
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

    for (int i = 0; i < length / 2; ++i) {
      tuplesMp[dataset[i].id] += dataset[i].value;
      tuplesMp1[dataset[i].id] += dataset[i].value;
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        each_value_total[j] += dataset[i].value.values[j];
      }
    }

    for (int i = length / 2; i < length; ++i) {
      tuplesMp[dataset[i].id] += dataset[i].value;
      tuplesMp2[dataset[i].id] += dataset[i].value;
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        each_value_total[j] += dataset[i].value.values[j];
      }
    }

    for (auto &[id, value] : tuplesMp) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        if (abs(tuplesMp1[id].values[j] - tuplesMp2[id].values[j]) >
            alpha * each_value_total[j]) {
          std::cout << alpha * each_value_total[j] << std::endl;
          std::cout << tuplesMp1[id].values[j] << " " << tuplesMp2[id].values[j]
                    << std::endl;
          heavychange[id].values[j] = 1;
        }
      }
    }
  }

  ~BenchMark_HeavyChange() {}

  void HHHyperBench(uint32_t MEMORY, uint32_t HASH_NUM = 2,
                    std::string BENCHNAME = "Ours") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start HyperBench!" << std::endl;

    OurHyper *sketch1 = new OurHyper(MEMORY, HASH_NUM);
    OurHyper *sketch2 = new OurHyper(MEMORY, HASH_NUM);

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (uint32_t i = 0; i < length / 2; ++i) {
      sketch1->Insert(dataset[i]);
    }
    for (uint32_t i = length / 2; i < length; ++i) {
      sketch2->Insert(dataset[i]);
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
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple1 = sketch1->AllQuery();
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple2 = sketch2->AllQuery();
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "Hyper Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

    std::set<TUPLES_ID> id_set;
    for (auto &[id, value] : estTuple1) {
      id_set.insert(id);
    }
    for (auto &[id, value] : estTuple2) {
      id_set.insert(id);
    }
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple;
    for (auto id : id_set) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        if (abs(estTuple1[id].values[j] - estTuple2[id].values[j]) >
            alpha * each_value_total[j]) {
          estTuple[id].values[j] = 1;
        }
      }
    }

    Compare_HeavyChange(estTuple, heavychange);

    delete sketch1;
    delete sketch2;
  }

  void HHSingleBench(uint32_t MEMORY, uint32_t HASH_NUM = 2,
                     std::string BENCHNAME = "Single") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start SingleBench!" << std::endl;

    OurSingle *sketch1[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch1[i] = new OurSingle(MEMORY / TUPLES_VALUES_ELEMENT_NUM);
    }

    OurSingle *sketch2[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch2[i] = new OurSingle(MEMORY / TUPLES_VALUES_ELEMENT_NUM);
    }

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < length / 2; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        Single tmp;
        tmp.id = dataset[i].id;
        tmp.value = dataset[i].value.values[j];
        sketch1[j]->Insert(tmp);
      }
    }
    for (int i = length / 2; i < length / 2; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        Single tmp;
        tmp.id = dataset[i].id;
        tmp.value = dataset[i].value.values[j];
        sketch2[j]->Insert(tmp);
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
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple1;
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple1(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple1[i] = sketch1[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple1[i].begin(); it != tmpTuple1[i].end(); ++it) {
        estTuple1[it->first].values[i] = it->second;
      }
    }
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple2;
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple2(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple2[i] = sketch2[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple2[i].begin(); it != tmpTuple2[i].end(); ++it) {
        estTuple2[it->first].values[i] = it->second;
      }
    }
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "Single Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

    std::set<TUPLES_ID> id_set;
    for (auto &[id, value] : estTuple1) {
      id_set.insert(id);
    }
    for (auto &[id, value] : estTuple2) {
      id_set.insert(id);
    }
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple;
    for (auto id : id_set) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        if (abs(estTuple1[id].values[j] - estTuple2[id].values[j]) >
            alpha * each_value_total[j]) {
          estTuple[id].values[j] = 1;
        }
      }
    }

    Compare_HeavyChange(estTuple, heavychange);

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      delete sketch1[i];
      delete sketch2[i];
    }
  }

  void HHUSSBench(uint32_t MEMORY, uint32_t HASH_NUM = 2,
                  std::string BENCHNAME = "USS") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start USSBench!" << std::endl;

    OurUSS *sketch1[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch1[i] = new OurUSS(MEMORY / TUPLES_VALUES_ELEMENT_NUM);
    }

    OurUSS *sketch2[TUPLES_VALUES_ELEMENT_NUM];

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      sketch2[i] = new OurUSS(MEMORY / TUPLES_VALUES_ELEMENT_NUM);
    }

    // insert
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < length / 2; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        sketch1[j]->Insert(dataset[i].id, dataset[i].value.values[j]);
      }
    }
    for (int i = length / 2; i < length; ++i) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; j++) {
        sketch2[j]->Insert(dataset[i].id, dataset[i].value.values[j]);
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
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple1;
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple1(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple1[i] = sketch1[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple1[i].begin(); it != tmpTuple1[i].end(); ++it) {
        estTuple1[it->first].values[i] = it->second;
      }
    }
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple2;
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple2(
        TUPLES_VALUES_ELEMENT_NUM);
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      tmpTuple2[i] = sketch2[i]->AllQuery();
    }
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      for (auto it = tmpTuple2[i].begin(); it != tmpTuple2[i].end(); ++it) {
        estTuple2[it->first].values[i] = it->second;
      }
    }
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tm_query =
        end_query - start_query;
    std::cout << "USS Query Time: " << tm_query.count() / 1000 << "s"
              << std::endl;

    std::set<TUPLES_ID> id_set;
    for (auto &[id, value] : estTuple1) {
      id_set.insert(id);
    }
    for (auto &[id, value] : estTuple2) {
      id_set.insert(id);
    }
    std::unordered_map<TUPLES_ID, TUPLES_VALUE> estTuple;
    for (auto id : id_set) {
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        if (abs(estTuple1[id].values[j] - estTuple2[id].values[j]) >
            alpha * each_value_total[j]) {
          estTuple[id].values[j] = 1;
        }
      }
    }

    Compare_HeavyChange(estTuple, heavychange);

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      delete sketch1[i];
      delete sketch2[i];
    }
  }

 private:
  std::string filename;
  std::string benchname;

  std::vector<TUPLES> dataset;

  double alpha;
  uint32_t length;
  uint32_t memory;
  double throughput;

  std::unordered_map<TUPLES_ID, TUPLES_VALUE> tuplesMp;
  std::unordered_map<TUPLES_ID, TUPLES_VALUE> tuplesMp1;
  std::unordered_map<TUPLES_ID, TUPLES_VALUE> tuplesMp2;
  std::unordered_map<TUPLES_ID, TUPLES_VALUE> heavychange;
  double each_value_total[TUPLES_VALUES_ELEMENT_NUM] = {0};

  template <class T>
  // mp : estimate
  // record : real
  void Compare_HeavyChange(T mp, T record) {
    std::set<TUPLES_ID> id_set;
    for (auto &[id, value] : mp) {
      id_set.insert(id);
    }
    for (auto &[id, value] : record) {
      id_set.insert(id);
    }
    double precision_sum = 0.0, recall_sum = 0.0;
    for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
      double bothHH = 0, realHH = 0, estHH = 0;
      for (auto id : id_set) {
        bool real = (record[id].values[j] == 1);
        bool est = (mp[id].values[j] == 1);
        realHH += real;
        estHH += est;
        if (real && est) {
          ++bothHH;
        }
      }
      std::cout << bothHH << " " << realHH << " " << estHH << std::endl;
      precision_sum += bothHH / estHH;
      recall_sum += bothHH / realHH;
    }
    double NUM = TUPLES_VALUES_ELEMENT_NUM;
    double F1 = (2 * precision_sum / NUM * recall_sum / NUM) /
                (precision_sum / NUM + recall_sum / NUM);
    std::ofstream ofs;
    ofs.open(filename + ".txt", std::ios::app);
    ofs << benchname << " " << memory << " " << precision_sum / NUM << " "
        << recall_sum / NUM << " " << F1 << std::endl;
    ofs.close();
  }
};

#endif