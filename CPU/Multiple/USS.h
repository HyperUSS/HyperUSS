#ifndef USS_H
#define USS_H

#include <cstring>
#include <random>
#include <unordered_map>

#include "StreamSummary.h"
#include "Util.h"

class OurUSS {
public:
  typedef std::unordered_map<TUPLES_ID, double> HashMap;

  OurUSS(uint32_t _MEMORY, std::string _name = "OurUSS") {
    this->name = _name;

    summary = new StreamSummary(summary->Memory2Size(_MEMORY));
  }

  ~OurUSS() { delete summary; }

  void Insert(const TUPLES_ID &item, double value) {
    if (summary->mp->Lookup(item)) {
      summary->Add_Data(item, value);
    } else {
      if (summary->isFull()) {
        double min_value = summary->getMin();
        double prob = value / (min_value + value);
        static std::mt19937 e2(rd());
        std::uniform_real_distribution<> dist(0, 1);
        if (dist(e2) < prob) {
          summary->SS_Replace(item, value);
        } else {
          summary->Add_Min(value);
        }
      } else {
        summary->New_Data(item, value);
      }
    }
  }

  HashMap AllQuery() { return summary->AllQuery(); }

private:
  std::string name;
  StreamSummary *summary;
};

#endif
