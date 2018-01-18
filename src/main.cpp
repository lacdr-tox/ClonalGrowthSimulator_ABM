#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>  // for high_resolution_clock
#include <ctime>
#include <map>
#include <random>

#include "SimulationBuilder.h"
#include "GrowthModel.h"
#include "Experiment.h"


int main(int argc, char *argv[]) {
  auto start = std::chrono::steady_clock::now();
  bool verbose = true;
  // find parameter file
  if (argc < 2) {
    std::cout << "Settings file missing\n";
    exit(EXIT_FAILURE);
  }
  if (argc == 3) {
    if ((std::string(argv[1])).compare("-q")) { verbose = false; }
    else if ((std::string(argv[1])).compare("-v")) { verbose = true; }

  }
  std::string fn = argv[argc - 1];
  SimulationBuilder builder = SimulationBuilder(fn,verbose);

  // Run simulation
  std::cout << "===== SIMULATION =====\n";
  builder.experiment->Run();
  std::cout << "===== POST-PROCESSING =====\n";
  if (builder.experiment->gzip)
    builder.experiment->GzipOutput();

  // Get simulation time
  if (verbose) {
    auto end = std::chrono::steady_clock::now();
    int tsec = (int) std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    int tmin = (int) std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
    int thour = (int) std::chrono::duration_cast<std::chrono::hours>(end - start).count();
    std::cout << "Simulation time: ";
    if (thour > 0)
      std::cout << thour << " hours ";
    if ((tmin > 0) || (thour > 0))
      std::cout << tmin - (thour * 60) << " minutes ";
    std::cout << tsec - (tmin * 60) << " seconds\n";
  }

  return 0;
}
