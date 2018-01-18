#ifndef ABM_SIMULATIONBUILDER_H
#define ABM_SIMULATIONBUILDER_H

#include <string>
#include "tinyxml2.h"
#include "Experiment.h"
#include "GrowthModel.h"

//! Simulation builder
/*!
 * \class SimulationBuilder
 *
 * The builder parses the xml, creates the model, solver and experiment
 * instances, and performs all actions needed before the experiment can run.
 */
class SimulationBuilder {
  public:
    GrowthModel *model;
    Experiment *experiment;

    SimulationBuilder();
    SimulationBuilder(std::string fn, bool verbose);

 private:
  void ToXML(std::string fn);
  std::map<std::string, std::string> ParsePars(const tinyxml2::XMLNode *node);
  void BuildModel(const tinyxml2::XMLNode *node);
  void BuildExperiment(const tinyxml2::XMLNode *node);
  void ParseFile(std::string fn);

};


#endif //ABM_SIMULATIONBUILDER_H
