#ifndef ABM_EXPERIMENTCONSTANTGROWTH_H
#define ABM_EXPERIMENTCONSTANTGROWTH_H

#include "Experiment.h"

class ExperimentConstantGrowth : public Experiment {

 public:
  //! Create an constant growth Experiment
  ExperimentConstantGrowth();

  //! Create an constant growth Experiment
  /*!
   *
   * \param pars map with parameters as read from the xml
   *
   */
  ExperimentConstantGrowth(std::map<std::string, std::string> pars);

  virtual ~ExperimentConstantGrowth();


  //! Run iterated growth experiment
  virtual void Run();

  //! Create map of simulation settings that can be used to generate the `Experiment` xml context
  virtual std::map<std::string, std::string> GetPars();

  //! simulation time
  double max_time;

  //! frequency for writing output
  int save_freq;

 private:
  double save_interval;
  bool save_final_cells,save_init_cells;

};


#endif //ABM_EXPERIMENTCONSTANTGROWTH_H
