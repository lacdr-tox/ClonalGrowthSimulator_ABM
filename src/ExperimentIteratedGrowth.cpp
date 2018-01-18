#include "ExperimentIteratedGrowth.h"
#include <iostream>
#include <boost/random/discrete_distribution.hpp>
#include <algorithm>
#include <string>
#include <iomanip>      // std::setprecision

ExperimentIteratedGrowth::ExperimentIteratedGrowth() : Experiment() {
  number_of_passages = 1;
  number_of_passing_cells = 1000;
  max_pop_size = 2000;
  max_pass_time = 1e9;
  grow_to_pop_size = true;
  stop_sim_if_n_pass_not_reached = false;
  biased_passage = false;
  type = "IteratedGrowth";
}

ExperimentIteratedGrowth::ExperimentIteratedGrowth(std::map<std::string, std::string> pars)
    : ExperimentIteratedGrowth() {
  SetPars(pars);
  for (auto par : pars) {
    if (par.first.compare("NumberOfCellsToKeep") == 0)
      number_of_passing_cells = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("CriticalPopulationSize") == 0)
      max_pop_size = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("NumberOfPassages") == 0)
      number_of_passages = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("MaxPassTime") == 0) {
      max_pass_time = std::stod(par.second);
      grow_to_pop_size = true;
    } else if (par.first.compare("StopSimIfNPassNotReached") == 0)
      stop_sim_if_n_pass_not_reached = true;
    else if (par.first.compare("BiasedPassage") == 0)
      biased_passage = true;
  }
  if (ncells_init == 0)
    ncells_init = number_of_passing_cells;
}

std::map<std::string, std::string> ExperimentIteratedGrowth::GetPars() {
  std::map<std::string, std::string> pars = GetBasePars();
  pars["NumberOfCellsToKeep"] = std::to_string(number_of_passing_cells);
  pars["NumberOfPassages"] = std::to_string(number_of_passages);
  pars["CriticalPopulationSize"] = std::to_string(max_pop_size);
  pars["MaxPassTime"] = std::to_string(max_pass_time);
  if (stop_sim_if_n_pass_not_reached) { pars["StopSimIfNPassNotReached"] = ""; }
  if (biased_passage) { pars["BiasedPassage"] = ""; }
  return pars;
};


void ExperimentIteratedGrowth::Run() {
  std::cout << "Run Iterated Growth\n";
  current_passage = 0;
  current_step = 0;
  double t = 0;
  SaveState(0, "clones_init");
  SavePopulationStats(0);
  if (track_divtime) { SaveDivRates(0, "divrates_init"); }
  if (track_divtime_bc) { SaveDivRatesBC(0, "divratesBC_init"); }
  while (current_passage < number_of_passages) {
    current_passage += 1;
    current_step += 1;
    std::cout << "Passage " << current_passage << std::endl;
    if (grow_to_pop_size) {
      while (model->nof_cells < max_pop_size) {
        t = model->Step();
        if (t >= current_passage * max_pass_time) {
          std::cout << "stop at t = " << t << std::endl;
          break;
        }
      }
    }
    SaveState(current_passage, "clones_before_passage");
    SavePopulationStats(t);
    if (model->nof_cells > number_of_passing_cells)
      Passage(model->gen);
    if (track_cells){SaveCellInfo(std::to_string(current_passage));}
    if (track_divtime) { SaveDivRates(current_passage); }
    if (track_divtime_bc) { SaveDivRatesBC(current_passage); }
    SaveState(current_passage, "clones_after_passage");
  }
}


void ExperimentIteratedGrowth::Passage(base_generator_type &gen) {
  if (biased_passage)
    PassageBiased(gen);
  else
    PassageUniform(gen);
};

void ExperimentIteratedGrowth::PassageUniform(base_generator_type &gen) {
  std::vector<int> idx;
  for (auto c : model->cells) { idx.push_back(c.first); }
  std::shuffle(idx.begin(), idx.end(), model->gen);
  std::vector<int> keep(idx.begin(), idx.begin() + number_of_passing_cells);
  model->UpdateQueues(keep);
};

void ExperimentIteratedGrowth::PassageBiased(base_generator_type &gen) {
  std::vector<double> probs;
  double num = (double) ((model->barcodes.size() - 1) * model->nof_cells);
  int cnt;
  for (auto c : model->cells) {
    cnt = model->barcodes[c.second.uid];
    probs.push_back((model->nof_cells - cnt) / (num * cnt));
  }
  boost::random::discrete_distribution<> dist(probs);
  std::set<int> keep;
  while ((int) keep.size() < number_of_passing_cells) {
    keep.insert(model->cells[dist(gen)].uid);
  }
  std::vector<int> vkeep(keep.begin(), keep.end());
  model->UpdateQueues(vkeep);
};
