#ifndef ABM_EXPERIMENTITERATEDGROWTH_H
#define ABM_EXPERIMENTITERATEDGROWTH_H

#include "Experiment.h"
#include <string>
#include <map>

class ExperimentIteratedGrowth : public Experiment  {
 public:
  ExperimentIteratedGrowth();
  ExperimentIteratedGrowth(std::map<std::string, std::string> pars);
  virtual ~ExperimentIteratedGrowth(){};

  //! Run iterated growth experiment
  virtual void Run();

  //! Create map of simulation settings that can be used to generate the `Experiment` xml context
  virtual std::map<std::string, std::string> GetPars();

  void Passage(base_generator_type &gen);
  void PassageUniform(base_generator_type &gen);
  void PassageBiased(base_generator_type &gen);

  //! number of passages performed during a simulation
  unsigned int number_of_passages;
  //! number of cells that is passed to the next generation
  int number_of_passing_cells;
  //! population size at which growth is stopped
  int max_pop_size;
  //! time (in days) at which growth is stopped
  double max_pass_time;
  //! stop simulating if the population size after growth < number of passing cells
  bool stop_sim_if_n_pass_not_reached;
  //! passage is biased by clone size
  bool biased_passage;
  unsigned int current_passage;
  bool grow_to_pop_size;
};


#endif //ABM_EXPERIMENTITERATEDGROWTH_H
