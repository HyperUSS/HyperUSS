#ifndef STREAMSUMMARY_H
#define STREAMSUMMARY_H

#include "CuckooMap.h"

template <typename ID_TYPE>
class Node {
 public:
  ID_TYPE ID;
  Node *prev;
  Node *next;

  Node(ID_TYPE _ID) : ID(_ID), prev(nullptr), next(nullptr) {}

  void Delete() { Connect(prev, next); }

  void Connect(Node *prev, Node *next) {
    if (prev) prev->next = next;
    if (next) next->prev = prev;
  }
};

class StreamSummary {
 public:
  typedef std::unordered_map<TUPLES_ID, double> HashMap;

  class DataNode;
  class CountNode;

  class DataNode : public Node<TUPLES_ID> {
   public:
    CountNode *pCount;
    DataNode(TUPLES_ID _ID) : Node<TUPLES_ID>(_ID), pCount(nullptr) {}
  };

  class CountNode : public Node<double> {
   public:
    DataNode *pData;
    CountNode(double _ID = 0) : Node<double>(_ID), pData(nullptr) {}
  };

  typedef CuckooMap<TUPLES_ID, DataNode *> Cuckoo;

  StreamSummary(uint32_t _SIZE) {
    SIZE = _SIZE;
    mp = new Cuckoo(SIZE);
    min = new CountNode();
  }

  ~StreamSummary() {
    delete mp;
    CountNode *pCount = min;
    while (pCount) {
      DataNode *pData = pCount->pData;
      while (pData) {
        DataNode *nextData = (DataNode *)pData->next;
        delete pData;
        pData = nextData;
      }

      CountNode *nextCount = (CountNode *)pCount->next;
      delete pCount;
      pCount = nextCount;
    }
  }

  static uint32_t Size2Memory(uint32_t size) {
    return size * ((sizeof(TUPLES_ID) + sizeof(DataNode *)) / LOAD +
                   sizeof(TUPLES_ID) + sizeof(double) + 4 * sizeof(void *));
  }

  static uint32_t Memory2Size(uint32_t memory) {
    return memory / ((sizeof(TUPLES_ID) + sizeof(DataNode *)) / LOAD +
                     sizeof(TUPLES_ID) + sizeof(double) + 4 * sizeof(void *));
  }

  uint32_t SIZE;
  Cuckoo *mp;
  CountNode *min;

  inline double getMin() { return min->next->ID; }

  double Query(const TUPLES_ID &item) {
    return mp->Lookup(item) ? (*mp)[item]->pCount->ID : 0;
  }

  inline bool isFull() { return mp->size() >= SIZE; }

  HashMap AllQuery() {
    HashMap ret;
    CountNode *pCount = min;
    while (pCount) {
      DataNode *pData = pCount->pData;
      while (pData) {
        ret[pData->ID] = pCount->ID;
        pData = (DataNode *)pData->next;
      }
      pCount = (CountNode *)pCount->next;
    }
    return ret;
  }

  inline void New_Data(const TUPLES_ID &data, double value) {
    DataNode *pData = new DataNode(data);
    Add_Count(min, pData, value);
    mp->Insert(data, pData);
  }

  inline void Add_Min(double value) {
    Add_Data(((CountNode *)min->next)->pData->ID, value);
  }

  void Add_Data(const TUPLES_ID &data, double value) {
    DataNode *pData = (*mp)[data];
    CountNode *pCount = pData->pCount;

    bool del = false;
    pData->Delete();
    if (pCount->pData == pData) {
      pCount->pData = (DataNode *)pData->next;
      del = !pData->next;
    }

    Add_Count(pCount, pData, value);

    if (del) {
      pCount->Delete();
      delete pCount;
    }
  }

  void Add_Count(CountNode *pCount, DataNode *pData, double insert_count) {
    CountNode *prev;
    CountNode *p;
    for (prev = pCount, p = (CountNode *)pCount->next; p;
         prev = p, p = (CountNode *)p->next) {
      if (p->ID - pCount->ID > insert_count) {
        CountNode *add = new CountNode(pCount->ID + insert_count);
        pCount->Connect(add, p);
        pCount->Connect(prev, add);
        p = add;
        break;
      } else if (p->ID == pCount->ID + insert_count)
        break;
    }
    if (!p) {
      pCount->Connect(prev, new CountNode(pCount->ID + insert_count));
      p = (CountNode *)prev->next;
    }

    pData->prev = nullptr;
    pData->pCount = (CountNode *)p;
    // pdata insert into new queue's header
    pData->Connect(pData, pData->pCount->pData);
    pData->pCount->pData = pData;
  }

  void SS_Replace(const TUPLES_ID &data, double value) {
    CountNode *pCount = (CountNode *)min->next;
    DataNode *pData = new DataNode(data);

    mp->Insert(data, pData);
    Add_Count(pCount, pData, value);
    pData = pCount->pData;
    pCount->pData = (DataNode *)pData->next;
    pData->Delete();
    if (!pData->next) {
      pCount->Delete();
      delete pCount;
    }
    mp->Delete(pData->ID);
    delete pData;
  }
};

#endif
