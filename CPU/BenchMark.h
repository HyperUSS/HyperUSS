#ifndef HHBENCH_H
#define HHBENCH_H

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

class BenchMark {
 public:
  BenchMark(std::string PATH, std::string FILENAME, double SKEW_FACTOR = 1) {
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
        double tmp;
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
        double tmp;
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
      for (int j = 0; j < TUPLES_VALUES_ELEMENT_NUM; ++j) {
        each_value_total[j] += dataset[i].value.values[j];
      }
    }
    std::cout << "Number of different full keys: " << tuplesMp.size()
              << std::endl;
  }

  ~BenchMark() {}

  void HHHyperBench(uint32_t MEMORY, double alpha, uint32_t HASH_NUM = 2,
                    std::string BENCHNAME = "Ours") {
    benchname = BENCHNAME;

    memory = MEMORY;

    std::cout << "start HyperBench!" << std::endl;

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
    // full key
    Compare_FullKey(estTuple, tuplesMp, alpha);

    // partial key
    Compare_PartialKey(estTuple, tuplesMp, alpha);
#else
    // full key
    Compare_FullKey_SKEW(estTuple, tuplesMp, alpha);

    // partial key
    Compare_PartialKey_SKEW(estTuple, tuplesMp, alpha);
#endif

    delete sketch;
  }

  void HHSingleBench(uint32_t MEMORY, double alpha, uint32_t HASH_NUM = 2,
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
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple(
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

    // full key
    Compare_FullKey(estTuple, tuplesMp, alpha);

    // partial key
    Compare_PartialKey(estTuple, tuplesMp, alpha);

    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; i++) {
      delete sketch[i];
    }
  }

  void HHUSSBench(uint32_t MEMORY, double alpha,
                  std::string BENCHNAME = "USS") {
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
    std::vector<std::unordered_map<TUPLES_ID, double>> tmpTuple(
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

    // full key
    Compare_FullKey(estTuple, tuplesMp, alpha);

    // partial key
    Compare_PartialKey(estTuple, tuplesMp, alpha);

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

  std::unordered_map<TUPLES_ID, TUPLES_VALUE> tuplesMp;
  double each_value_total[TUPLES_VALUES_ELEMENT_NUM] = {0};

  template <class T>
  // mp : estimate
  // record : real
  void Compare_FullKey(T mp, T record, double alpha) {
    double precision_sum = 0.0;
    double recall_sum = 0.0;
    double aae_sum = 0.0;
    double are_sum = 0.0;
#ifndef AVERAGE
    double NUM = TUPLES_VALUES_ELEMENT_NUM;
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      double threshold = each_value_total[i] * alpha;
      double bothHH = 0, realHH = 0, estHH = 0;
      double aae = 0, are = 0;
      for (auto it = record.begin(); it != record.end(); ++it) {
        double realF = it->second.values[i];
        double estF = mp[it->first].values[i];
        bool real = (abs(realF) > threshold);
        bool est = (abs(estF) > threshold);
        realHH += real;
        estHH += est;
        if (real && est && realF != 0) {
          bothHH += 1;
          aae += abs(realF - estF);
          are += abs(realF - estF) / abs(realF);
        }
      }
      if (bothHH != 0) {
        precision_sum += bothHH / estHH;
        recall_sum += bothHH / realHH;
        aae_sum += aae / bothHH;
        are_sum += are / bothHH;
      } else if (estHH == 0 && realHH == 0) {
        precision_sum += 1.0;
        recall_sum += 1.0;
      }
    }
#else
    double NUM = TUPLES_VALUES_ELEMENT_NUM - 1;
    for (int i = 1; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      double threshold = abs(each_value_total[0] * alpha);
      double bothHH = 0, realHH = 0, estHH = 0;
      double aae = 0, are = 0;
      for (auto it = record.begin(); it != record.end(); ++it) {
        double real_cnt = it->second.values[0];
        double est_cnt = mp[it->first].values[0];
        double realF = it->second.values[i];
        double estF = mp[it->first].values[i];
        bool real = (real_cnt > threshold);
        bool est = (est_cnt > threshold);
        realHH += real;
        estHH += est;
        if (real && est && realF != 0) {
          bothHH += 1;
          realF /= real_cnt;
          estF /= est_cnt;
          aae += abs(realF - estF);
          are += abs(realF - estF) / abs(realF);
        }
      }
      if (bothHH != 0) {
        precision_sum += bothHH / estHH;
        recall_sum += bothHH / realHH;
        aae_sum += aae / bothHH;
        are_sum += are / bothHH;
      } else if (estHH == 0 && realHH == 0) {
        precision_sum += 1.0;
        recall_sum += 1.0;
      }
    }
#endif
    std::ofstream ofs;
    double F1 = (2 * precision_sum / NUM * recall_sum / NUM) /
                (precision_sum / NUM + recall_sum / NUM);
    ofs.open(filename + "_FullKey.txt", std::ios::app);
    ofs << benchname << " " << memory << " " << precision_sum / NUM << " "
        << recall_sum / NUM << " " << F1 << " " << aae_sum / NUM << " "
        << are_sum / NUM << " " << throughput << std::endl;
    ofs.close();
  }

  template <class T>
  // mp : estimate
  // record : real
  void Compare_PartialKey(T mp, T record, double alpha) {
    double precision_sum = 0.0;
    double recall_sum = 0.0;
    double aae_sum = 0.0;
    double are_sum = 0.0;
    std::vector<std::vector<uint32_t>> keylist;
    std::ifstream ifs("../PartialKey.txt");
    if (!ifs) {
      std::cout << "fail to open PartialKey" << std::endl;
      abort();
    }
    std::vector<uint32_t> tmplist;
    while (!ifs.eof()) {
      uint32_t tmp;
      ifs >> tmp;
      tmplist.push_back(tmp);
      if (ifs.peek() == '\n' || ifs.eof()) {
        keylist.push_back(tmplist);
        std::vector<uint32_t>().swap(tmplist);
      }
    }
    for (int k = 0; k < keylist.size(); ++k) {
      std::unordered_map<TUPLES_ID, TUPLES_VALUE> tmp_mp;
      std::unordered_map<TUPLES_ID, TUPLES_VALUE> tmp_record;
      for (auto &[x, y] : mp) {
        TUPLES_ID id = TUPLES_ID();
        for (auto p : keylist[k]) {
          id.key[p - 1] = x.key[p - 1];
        }
        tmp_mp[id] += y;
      }
      for (auto &[x, y] : record) {
        TUPLES_ID id;
        for (auto p : keylist[k]) {
          id.key[p - 1] = x.key[p - 1];
        }
        tmp_record[id] += y;
      }
#ifndef AVERAGE
      for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        double threshold = abs(each_value_total[i] * alpha);
        double bothHH = 0, realHH = 0, estHH = 0;
        double aae = 0, are = 0;
        for (auto it = tmp_record.begin(); it != tmp_record.end(); ++it) {
          double realF = it->second.values[i];
          double estF = tmp_mp[it->first].values[i];
          bool real = (abs(realF) > threshold);
          bool est = (abs(estF) > threshold);
          realHH += real;
          estHH += est;
          if (real && est && realF != 0) {
            bothHH += 1;
            aae += abs(realF - estF);
            are += abs(realF - estF) / abs(realF);
          }
        }
        if (bothHH != 0) {
          precision_sum += bothHH / estHH;
          recall_sum += bothHH / realHH;
          aae_sum += aae / bothHH;
          are_sum += are / bothHH;
        } else if (estHH == 0 && realHH == 0) {
          precision_sum += 1.0;
          recall_sum += 1.0;
        }
      }
#else
      for (int i = 1; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        double threshold = each_value_total[0] * alpha;
        double bothHH = 0, realHH = 0, estHH = 0;
        double aae = 0, are = 0;
        for (auto it = tmp_record.begin(); it != tmp_record.end(); ++it) {
          double real_cnt = it->second.values[0];
          double est_cnt = tmp_mp[it->first].values[0];
          double realF = it->second.values[i];
          double estF = tmp_mp[it->first].values[i];
          bool real = (real_cnt > threshold);
          bool est = (est_cnt > threshold);
          realHH += real;
          estHH += est;
          if (real && est && realF > 0) {
            bothHH += 1;
            realF /= real_cnt;
            estF /= est_cnt;
            aae += abs(realF - estF);
            are += abs(realF - estF) / abs(realF);
          }
        }
        precision_sum += bothHH / estHH;
        recall_sum += bothHH / realHH;
        aae_sum += aae / bothHH;
        are_sum += are / bothHH;
      }
#endif
    }
#ifndef AVERAGE
    double NUM = keylist.size() * TUPLES_VALUES_ELEMENT_NUM;
#else
    double NUM = keylist.size() * (TUPLES_VALUES_ELEMENT_NUM - 1);
#endif
    std::ofstream ofs;
    double F1 = (2 * precision_sum / NUM * recall_sum / NUM) /
                (precision_sum / NUM + recall_sum / NUM);
    ofs.open(filename + "_PartialKey.txt", std::ios::app);
    ofs << benchname << " " << memory << " " << precision_sum / NUM << " "
        << recall_sum / NUM << " " << F1 << " " << aae_sum / NUM << " "
        << are_sum / NUM << " " << throughput << std::endl;
    ofs.close();
  }

  template <class T>
  // mp : estimate
  // record : real
  void Compare_FullKey_SKEW(T mp, T record, double alpha) {
    double are_sum = 0.0;
    double min_are = std::numeric_limits<double>::infinity();
    double max_are = 0;
    for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
      double threshold = each_value_total[i] * alpha;
      double bothHH = 0, realHH = 0, estHH = 0;
      double aae = 0, are = 0;
      for (auto it = record.begin(); it != record.end(); ++it) {
        double realF = it->second.values[i];
        double estF = mp[it->first].values[i];
        bool real = (realF > threshold);
        bool est = (estF > threshold);
        realHH += real;
        estHH += est;
        if (real && est) {
          bothHH += 1;
          aae += abs(realF - estF);
          are += abs(realF - estF) / abs(realF);
        }
      }
      are_sum += are / bothHH;
      min_are = are / bothHH < min_are ? are / bothHH : min_are;
      max_are = are / bothHH > max_are ? are / bothHH : max_are;
    }
    std::ofstream ofs;
    ofs.open(filename + "_FullKey.txt", std::ios::app);
    ofs << memory << " " << skew_factor << " " << min_are << " " << max_are
        << " " << are_sum / TUPLES_VALUES_ELEMENT_NUM << " " << std::endl;
    ofs.close();
  }

  template <class T>
  // mp : estimate
  // record : real
  void Compare_PartialKey_SKEW(T mp, T record, double alpha) {
    std::vector<std::vector<uint32_t>> keylist;
    std::ifstream ifs("../PartialKey.txt");
    if (!ifs) {
      std::cout << "fail to open PartialKey" << std::endl;
      abort();
    }
    std::vector<uint32_t> tmplist;
    while (!ifs.eof()) {
      uint32_t tmp;
      ifs >> tmp;
      tmplist.push_back(tmp);
      if (ifs.peek() == '\n' || ifs.eof()) {
        keylist.push_back(tmplist);
        std::vector<uint32_t>().swap(tmplist);
      }
    }
    std::vector<double> are_vec(TUPLES_VALUES_ELEMENT_NUM, 0);
    for (int k = 0; k < keylist.size(); ++k) {
      std::unordered_map<TUPLES_ID, TUPLES_VALUE> tmp_mp;
      std::unordered_map<TUPLES_ID, TUPLES_VALUE> tmp_record;
      for (auto &[x, y] : mp) {
        TUPLES_ID id = TUPLES_ID();
        for (auto p : keylist[k]) {
          id.key[p - 1] = x.key[p - 1];
        }
        tmp_mp[id] += y;
      }
      for (auto &[x, y] : record) {
        TUPLES_ID id;
        for (auto p : keylist[k]) {
          id.key[p - 1] = x.key[p - 1];
        }
        tmp_record[id] += y;
      }
      for (int i = 0; i < TUPLES_VALUES_ELEMENT_NUM; ++i) {
        double threshold = each_value_total[i] * alpha;
        double bothHH = 0, realHH = 0, estHH = 0;
        double aae = 0, are = 0;
        for (auto it = tmp_record.begin(); it != tmp_record.end(); ++it) {
          double realF = it->second.values[i];
          double estF = tmp_mp[it->first].values[i];
          bool real = (realF > threshold);
          bool est = (estF > threshold);
          realHH += real;
          estHH += est;
          if (real && est) {
            bothHH += 1;
            aae += abs(realF - estF);
            are += abs(realF - estF) / abs(realF);
          }
        }
        are_vec[i] += are / bothHH;
      }
    }
    std::ofstream ofs;
    ofs.open(filename + "_PartialKey.txt", std::ios::app);
    ofs << memory << " " << skew_factor << " "
        << *std::min_element(are_vec.begin(), are_vec.end()) / keylist.size()
        << " "
        << *std::max_element(are_vec.begin(), are_vec.end()) / keylist.size()
        << " "
        << std::accumulate(are_vec.begin(), are_vec.end(), 0.0) /
               (keylist.size() * are_vec.size())
        << std::endl;
    ofs.close();
  }
};

#endif
