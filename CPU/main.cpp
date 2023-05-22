#include "BenchMark.h"
#include "BenchMark_HeavyChange.h"
#include "BenchMark_Synthetic.h"

int main(int argc, char *argv[]) {
  // example
  for (uint32_t i = 1; i < argc; ++i) {
    std::cout << "example dataset: " << argv[i] << std::endl;
    BenchMark_Synthetic dataset(argv[i], "example");
    dataset.HHHyperBench(1000 * 500);
    dataset.HHSingleBench(1000 * 500);
    dataset.HHUSSBench(1000 * 500);
  }

  return 0;
}

// experiments:

// exp1a
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp1a dataset: " << argv[i] << std::endl;
//     BenchMark_Synthetic dataset(argv[i], "exp1a");
//     for (int j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(1000 * j);
//       dataset.HHSingleBench(1000 * j);
//       dataset.HHUSSBench(1000 * j);
//     }
//   }

// exp1b
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp1b dataset: " << argv[i] << std::endl;
//     BenchMark_Synthetic dataset(argv[i], "exp1b");
//     for (int j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(1000 * j);
//       dataset.HHSingleBench(1000 * j);
//       dataset.HHUSSBench(1000 * j);
//     }
//   }

// exp1c
//   int i = 5;
//   std::string path =
//       "../data/Synthetic/Hyper_zipf_1.500_5K" + std::to_string(i) +
//       "V.txt";
//   std::cout << "exp1c dataset: " << path << std::endl;
//   BenchMark_Synthetic dataset(path, "exp1c");
//   dataset.HHHyperBench(1000 * 500);
//   dataset.HHSingleBench(1000 * 500);
//   dataset.HHUSSBench(1000 * 500);

// exp1d
//   int i = 5;
//   std::string path =
//       "../data/Synthetic/Hyper_zipf_1.500_" + std::to_string(i) +
//       "K5V.txt";
//   std::cout << "exp1d dataset: " << path << std::endl;
//   BenchMark_Synthetic dataset(path, "exp1d");
//   dataset.HHHyperBench(1000 * 500);
//   dataset.HHSingleBench(1000 * 500);
//   dataset.HHUSSBench(1000 * 500);

// exp1e
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp1e dataset: " << argv[i] << std::endl;
//     BenchMark_Synthetic dataset(argv[i], "exp1e");
//     for (int j = 1; j <= 5; j += 1) {
//       dataset.HHHyperBench(1000 * 500, j);
//     }
//   }

// exp1f
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp1f dataset: " << argv[i] << std::endl;
//     BenchMark_Synthetic dataset(argv[i], "exp1f");
//     for (int j = 800; j <= 1200; j += 100) {
//       dataset.HHHyperBench(1000 * j);
//       dataset.HHSingleBench(1000 * j);
//       dataset.HHUSSBench(1000 * j);
//     }
//   }

// exp1g
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp1g dataset: " << argv[i] << std::endl;
//     BenchMark_Synthetic dataset(argv[i], "exp1g");
//     for (int j = 800; j <= 1200; j += 100) {
//       dataset.HHHyperBench(1000 * j);
//       dataset.HHSingleBench(1000 * j);
//       dataset.HHUSSBench(1000 * j);
//     }
//   }

// exp2a
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp2a dataset: " << argv[i] << std::endl;
//     BenchMark dataset(argv[i], "exp2a");
//     for (int j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(1000 * j, 0.000001);
//       dataset.HHSingleBench(1000 * j, 0.000001);
//       dataset.HHUSSBench(1000 * j, 0.000001);
//     }
//   }

// exp2b
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp2b dataset: " << argv[i] << std::endl;
//     BenchMark dataset(argv[i], "exp2b");
//     for (int j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(1000 * j, 0.000001);
//       dataset.HHSingleBench(1000 * j, 0.000001);
//       dataset.HHUSSBench(1000 * j, 0.000001);
//     }
//   }

// exp2c skew_factor = 1e5
// for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "exp2c dataset: " << argv[i] << std::endl;
//     BenchMark dataset(argv[i], "exp2c_1e5", 1e5);
//     for (int k = 0; k < 1000; ++k) {
//       dataset.HHHyperBench(1000 * 500, 0.000001);
//     }
//   }

// revision: NBA
// for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "revision dataset: " << argv[i] << std::endl;
//     BenchMark dataset(argv[i], "revision_NBA");
//     for (uint32_t j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(100000 * j, 0.000001);
//       dataset.HHSingleBench(100000 * j, 0.000001);
//       dataset.HHUSSBench(100000 * j, 0.000001);
//     }
//   }

// revision: CAIDA
// for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "revision dataset: " << argv[i] << std::endl;
//     BenchMark dataset(argv[i], "revision_CAIDA_FAIR");
//     for (uint32_t k = 0; k < 100; ++k) {
//       dataset.HHHyperBench(100000 * 500, 0.000001);
//     }
//   }

// revision: valuedistribution
// for (uint32_t i = 1; i < argc; ++i) {
//   std::cout << "revision dataset: " << argv[i] << std::endl;
//   BenchMark dataset(argv[i], "revision_isclick");
//   for (uint32_t j = 300; j <= 700; j += 100) {
//     dataset.HHHyperBench(10000 * j, 0);
//     dataset.HHSingleBench(10000 * j, 0);
//     dataset.HHUSSBench(10000 * j, 0);
//   }
// }

// revision: heavychange
//   for (uint32_t i = 1; i < argc; ++i) {
//     std::cout << "revision dataset: " << argv[i] << std::endl;
//     BenchMark_HeavyChange dataset(argv[i], "revision_heavychange");
//     for (uint32_t j = 300; j <= 700; j += 100) {
//       dataset.HHHyperBench(1000 * j);
//       dataset.HHSingleBench(1000 * j);
//       dataset.HHUSSBench(1000 * j);
//     }
//   }
