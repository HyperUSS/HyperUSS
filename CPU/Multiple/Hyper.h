#ifndef HYPER_H
#define HYPER_H

#include <cstring>
#include <random>
#include <unordered_map>

#include "Util.h"

typedef std::unordered_map<TUPLES_ID, TUPLES_VALUE> HashMap;

class OurHyper {
 public:
  OurHyper(uint32_t _MEMORY, uint32_t _HASH_NUM = 2,
           std::string _name = "OurHyper") {
    this->name = _name;

    HASH_NUM = _HASH_NUM;

    LENGTH = _MEMORY / _HASH_NUM / sizeof(TUPLES);

    counter = new TUPLES *[HASH_NUM];

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      counter[i] = new TUPLES[LENGTH];
      memset(counter[i], 0, sizeof(TUPLES) * LENGTH);
    }
  }

  ~OurHyper() {
    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      delete[] counter[i];
    }
    delete[] counter;
  }

  void Insert(const TUPLES &item) {
    double minimum = std::numeric_limits<double>::infinity();
    uint32_t minPos, minHash;

#ifdef FAIRNESS
    S += item.value;
#endif

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      uint32_t position = hash(item.id, i) % LENGTH;
      if (counter[i][position].id == item.id) {
        counter[i][position].value += item.value;
        return;
      }
#ifdef FAIRNESS
      if (counter[i][position].value.normalize(S).sum_squares() < minimum) {
        minPos = position;
        minHash = i;
        minimum = counter[i][position].value.normalize(S).sum_squares();
      }
#else
      if (counter[i][position].value.sum_squares() < minimum) {
        minPos = position;
        minHash = i;
        minimum = counter[i][position].value.sum_squares();
      }
#endif
    }
#ifdef FAIRNESS
    double sqrt1 = sqrt(item.value.normalize(S).sum_squares());
    double sqrt2 =
        sqrt(counter[minHash][minPos].value.normalize(S).sum_squares());
#else
    double sqrt1 = sqrt(item.value.sum_squares());
    double sqrt2 = sqrt(counter[minHash][minPos].value.sum_squares());
#endif
    if (sqrt1 + sqrt2 != 0) {
      double prob = sqrt1 / (sqrt1 + sqrt2);
      static std::mt19937 e2(rd());
      std::uniform_real_distribution<> dist(0, 1);

      if (dist(e2) < prob) {
        counter[minHash][minPos] = item;
        counter[minHash][minPos].value.divide(prob);
      } else {
        counter[minHash][minPos].value.divide(1 - prob);
      }
    }
  }

  HashMap AllQuery() {
    HashMap ret;

    uint32_t empty = 0;

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      for (uint32_t j = 0; j < LENGTH; ++j) {
        if (counter[i][j].id == TUPLES_ID() && counter[i][j].value.empty()) {
          ++empty;
        } else {
          ret[counter[i][j].id] += counter[i][j].value;
        }
      }
    }

    return ret;
  }

 private:
  std::string name;
  uint32_t LENGTH;
  uint32_t HASH_NUM;

  TUPLES **counter;

#ifdef FAIRNESS
  TUPLES_VALUE S = TUPLES_VALUE();
#endif
};

#endif
