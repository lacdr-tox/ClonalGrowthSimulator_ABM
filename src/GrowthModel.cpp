#include <iostream>
#include <cmath>
#include "GrowthModel.h"
#include <algorithm>    // std::random_shuffle

GrowthModel::GrowthModel() {
  division_rate = 1;
  death_rate = 0;
  division_rate_sd = 0;
  mutation_sd = 0;
}

GrowthModel::GrowthModel(double _division_rate, double _death_rate, double _division_rate_sd, double _mutation_sd){
  division_rate = _division_rate;
  death_rate = _death_rate;
  division_rate_sd = _division_rate_sd;
  mutation_sd = _mutation_sd;
  minimum_division_rate = 0;
  maximum_division_rate = -1;
  mutation_distr_gauss = std::normal_distribution<>(1,mutation_sd);
}

GrowthModel::GrowthModel(std::map<std::string, std::string> pars) : GrowthModel(){
  for (auto par : pars) {
    if (par.first.compare("DivisionRate") == 0)
      division_rate = std::stod(par.second);
    else if (par.first.compare("DeathRate") == 0)
      death_rate = std::stod(par.second);
    else if (par.first.compare("MutationSD") == 0)
      mutation_sd = std::stod(par.second);
    else if (par.first.compare("DivisionRateSD") == 0)
      division_rate_sd = std::stod(par.second);
    else if (par.first.compare("MinDivisionRate") == 0)
      minimum_division_rate = std::stod(par.second);
    else if (par.first.compare("MaxDivisionRate") == 0)
      maximum_division_rate = std::stod(par.second);
  }
  mutation_distr_gauss = std::normal_distribution<>(1,mutation_sd);
}

std::map<std::string, std::string> GrowthModel::GetPars() {
  std::map<std::string, std::string> pars = {{"DivisionRate", std::to_string(division_rate)},
                                             {"DeathRate", std::to_string(death_rate)},
                                             {"DivisionRateSD", std::to_string(division_rate_sd)},
                                             {"MutationSD", std::to_string(mutation_sd)},
                                             {"MinDivisionRate", std::to_string(minimum_division_rate)}};
  if (maximum_division_rate > 0)
    pars["MaxDivisionRate"] = maximum_division_rate;
  return pars;
}

void GrowthModel::PrintInfo() {
  std::map<std::string, std::string> pars = this->GetPars();
  for (auto par : pars){
    std::cout << par.first << " = " << par.second << std::endl;
  }
};

void GrowthModel::AddInitCell(int uid, int bc, double div_rate){
  cells[uid] = Cell(div_rate, death_rate, bc, uid);
  cells[uid].div_time = cells[uid].GetDivTime(gen);
  div_queue.insert(std::pair<double, int>(cells[uid].div_time,uid ));
  if (death_rate > 0){
    cells[uid].death_time = cells[uid].GetDeathTime(gen);
    death_queue.insert(std::pair<double, int>(cells[uid].death_time,uid ));
  }
  barcodes[bc]++;
}

void GrowthModel::AddCells(int ncells){
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  mutation_distr_gauss = std::normal_distribution<>(1,division_rate_sd);
  for (int uid = 0; uid < ncells; uid++){
    cells[uid] = Cell(division_rate, death_rate, -1, uid);
    if (division_rate_sd > 0){this->Mutate(cells[uid]);}
    cells[uid].div_time = cells[uid].GetDivTime(gen);
    div_queue.insert(std::pair<double, int>(cells[uid].div_time,uid ));
    if (death_rate > 0){
      cells[uid].death_time = cells[uid].GetDeathTime(gen);
      death_queue.insert(std::pair<double, int>(cells[uid].death_time,uid )); }
  }
  nof_cells = (int)cells.size();
  max_uid = nof_cells;
  mutation_distr_gauss = std::normal_distribution<>(1,mutation_sd);
}


void GrowthModel::AssignBarcodes(std::vector<int> barcode_distribution){
  barcodes.clear();
  std::vector<int> bc_list;
  for (unsigned int bc = 0; bc < barcode_distribution.size(); bc++) {
    barcodes[bc] = 0;
    for (int i = 0; i < barcode_distribution[bc]; i++) {bc_list.push_back(bc);}
  }
  std::vector<int> uid_list;
  for (auto item : cells){uid_list.push_back(item.first);}
  std::shuffle(uid_list.begin(), uid_list.end(), gen);
  for (unsigned int i = 0; i < bc_list.size(); i++){
    cells[uid_list[i]].barcode = bc_list[i];
    barcodes[bc_list[i]]++;
  }
}

double GrowthModel::Step(){
  Cell c;
  double t = 0;
  if ((death_queue.size() == 0) || (div_queue.begin()->first < death_queue.begin()->first)){
    t = div_queue.begin()->first;
    c = cells[div_queue.begin()->second];
    div_queue.erase(div_queue.begin());
    if (c.dead){
      this->Step();
    }
    else {
      this->Divide(c,t);
    }
  }
  else {
    t = death_queue.begin()->first;
    c = cells[death_queue.begin()->second];
    death_queue.erase(death_queue.begin());
    Death(c);
  }
  return t;
}

void GrowthModel::Death(Cell & c){
  c.dead = true;
  barcodes[c.barcode] -= 1;
  nof_cells -= 1;
  cells.erase(c.uid);
}

void GrowthModel::Divide(Cell & c, double t) {
  for (int i = 0; i < 2; i++){
    Cell daugther = Cell(c.division_rate, c.death_rate, c.barcode, max_uid+1);
    daugther.parent_uid = c.uid;
    daugther.birth_time = c.div_time;
    if (mutation_sd > 0){ this->Mutate(daugther);}
    cells[max_uid+1] = daugther;
    cells[max_uid+1].div_time = cells[max_uid+1].GetDivTime(gen)+t;
    div_queue.insert(std::pair<double,int>(cells[max_uid+1].div_time,daugther.uid));
    if (death_rate > 0){
      cells[max_uid+1].death_time = cells[max_uid+1].GetDeathTime(gen)+t;
      death_queue.insert(std::pair<double, int>(cells[max_uid+1].death_time,daugther.uid)); }
    max_uid += 1;
  }
  cells.erase(c.uid);
  nof_cells = (int)cells.size();
  barcodes[c.barcode] += 1;
}

void GrowthModel::Mutate(Cell &c) {
  double change = mutation_distr_gauss(gen);
  c.division_rate = change > minimum_division_rate ? change * c.division_rate : 0;
  if ((maximum_division_rate > 0) && (c.division_rate > maximum_division_rate))
    c.division_rate = maximum_division_rate;
}

void GrowthModel::UpdateQueues(std::vector<int> keep){
  div_queue.clear();
  death_queue.clear();
  std::map<int, Cell> new_cells;
  Cell c;
  for (auto bc : barcodes){ barcodes[bc.first] = 0;}
  for (auto uid : keep){
    c = cells[uid];
    barcodes[c.barcode] += 1;
    new_cells[uid] = c;
    div_queue.insert(std::pair<double, int>(c.div_time,uid));
    if (c.death_time > 0)
      death_queue.insert(std::pair<double, int>(c.death_time,uid));
  }
  cells = new_cells;
  nof_cells = (int)div_queue.size();
}
