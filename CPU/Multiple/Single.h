#ifndef SINGLE_H
#define SINGLE_H

#include <cstring>
#include <random>
#include <unordered_map>

#include "Util.h"

struct Single {
  TUPLES_ID id;
  double value;
};

typedef std::unordered_map<TUPLES_ID, double> HashMap_Single;

class OurSingle {
 public:
  OurSingle(uint32_t _MEMORY, uint32_t _HASH_NUM = 2,
            std::string _name = "OurSingle") {
    this->name = _name;

    HASH_NUM = _HASH_NUM;

    LENGTH = _MEMORY / _HASH_NUM / sizeof(Single);

    counter = new Single *[HASH_NUM];

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      counter[i] = new Single[LENGTH];
      memset(counter[i], 0, sizeof(Single) * LENGTH);
    }
  }

  ~OurSingle() {
    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      delete[] counter[i];
    }
    delete[] counter;
  }

  void Insert(const Single &item) {
    double minimum = std::numeric_limits<double>::infinity();
    uint32_t minPos, minHash;

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      uint32_t position = hash(item.id, i) % LENGTH;
      if (counter[i][position].id == item.id) {
        counter[i][position].value += item.value;
        return;
      }
      if (abs(counter[i][position].value) < minimum) {
        minPos = position;
        minHash = i;
        minimum = abs(counter[i][position].value);
      }
    }
    if (item.value == 0) return;
    counter[minHash][minPos].value += item.value;
    static std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    if (dist(e2) < abs(item.value) / abs(counter[minHash][minPos].value)) {
      counter[minHash][minPos].id = item.id;
    }
  }

  HashMap_Single AllQuery() {
    HashMap_Single ret;

    for (uint32_t i = 0; i < HASH_NUM; ++i) {
      for (uint32_t j = 0; j < LENGTH; ++j) {
        ret[counter[i][j].id] = counter[i][j].value;
      }
    }

    return ret;
  }

 private:
  std::string name;
  uint32_t LENGTH;
  uint32_t HASH_NUM;

  Single **counter;
};

#endif
