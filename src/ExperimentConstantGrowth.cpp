#include "ExperimentConstantGrowth.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <set>
#include <iostream>

ExperimentConstantGrowth::ExperimentConstantGrowth() : Experiment() {
  max_time = 1;
  save_interval = max_time;
  type = "ConstantGrowth";
  save_final_cells = false;
}

ExperimentConstantGrowth::ExperimentConstantGrowth(std::map<std::string, std::string> pars) :
    ExperimentConstantGrowth() {
  SetPars(pars);
  for (auto par : pars) {
    if (par.first.compare("SimulationTime") == 0)
      max_time = std::stod(par.second);
    else if (par.first.compare("SaveInterval") == 0)
      save_interval = std::stod(par.second);
    else if (par.first.compare("SaveFinalCells") == 0)
      save_final_cells = true;
    else if (par.first.compare("SaveInitCells") == 0)
      save_init_cells = true;
  }
}

std::map<std::string, std::string> ExperimentConstantGrowth::GetPars() {
  std::map<std::string, std::string> pars = GetBasePars();
  pars["SimulationTime"] = std::to_string(max_time);
  pars["SaveInterval"] = std::to_string(save_interval);
  if (save_final_cells){ pars["SaveFinalCells"] = "";}
  if (save_init_cells){ pars["SaveInitCells"] = "";}

  return pars;
}


ExperimentConstantGrowth::~ExperimentConstantGrowth() { }


void ExperimentConstantGrowth::Run() {
  std::cout << "Run constant growth for " << max_time << " days\n";
  current_time = 0;
  // next percentage of max_time to write output to stdout
  int p_next = 0;
  // next time to write data to file
  double s_next = save_interval;
  double last_save = 0;
  // write initial data to file
  SaveState(0);
  SavePopulationStats(0);
  if (save_init_cells){SaveCellInfo("init");}
  if ((track_divtime) || (model->division_rate_sd > 0)){ SaveDivRates(0,"divrates_init");}
  if (track_divtime_bc){ SaveDivRatesBC(0,"divratesBC_init");}
  // start growth simulation
  while (current_time < max_time) {
    current_time = model->Step();
    // write stats to stdout
    if ((100*current_time)/max_time > p_next){
      std::cout << p_next << "% - " << model->nof_cells << " cells\n";
      p_next += 10;
    }
    // write data to file
    if (current_time >= s_next){
      SaveState(current_time);
      SavePopulationStats(current_time);
      if (track_divtime){ SaveDivRates(current_time);}
      if (track_divtime_bc){ SaveDivRatesBC(current_time);}
      s_next += save_interval;
      last_save = current_time;
    }
  }
  // write state from the last step to file
  if (current_time != last_save) {
    SaveState(current_time);
    SavePopulationStats(current_time);
  }
  if (save_final_cells){SaveCellInfo("final");}

}
