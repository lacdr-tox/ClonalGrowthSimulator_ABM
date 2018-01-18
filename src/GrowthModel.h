#ifndef ABM_GROWTHMODEL_H
#define ABM_GROWTHMODEL_H

#include <map>
#include <string>
#include <map>
#include <vector>
#include <random>
#include "Cell.h"
typedef	 std::mt19937 base_generator_type;

class GrowthModel {
 public:
  GrowthModel();
  GrowthModel(double _division_rate, double _death_rate, double _division_rate_sd, double mutation_sd);
  GrowthModel(std::map<std::string, std::string> pars);
  void AddCells(int ncells);
  void AddInitCell(int uid, int bc, double div_rate);

  void AssignBarcodes(std::vector<int> barcode_distribution);
  double Step();
  base_generator_type gen;
  std::map<std::string, std::string> GetPars();
  void UpdateQueues(std::vector<int> keep);
  void PrintInfo();
  int nof_cells;
  int max_uid;
  std::map<int, int> barcodes;
  // cells contains all cells that are current present
  std::map<int, Cell> cells;
  int nclones0,real_nclones0;
  double division_rate_sd;
  double division_rate, death_rate;

 private:
  void Divide(Cell &c, double t);
  void Death(Cell &c);
  void Mutate(Cell &c);
  std::multimap<double, int> div_queue;
  std::multimap<double, int> death_queue;
  double mutation_sd;
  double minimum_division_rate,maximum_division_rate;
  std::normal_distribution<> mutation_distr_gauss;
};


#endif //ABM_GROWTHMODEL_H
