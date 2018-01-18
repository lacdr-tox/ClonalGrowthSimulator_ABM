#include "SimulationBuilder.h"
#include "ExperimentConstantGrowth.h"
#include "ExperimentIteratedGrowth.h"
#include <time.h>
#include <iostream>
using namespace tinyxml2;


SimulationBuilder::SimulationBuilder() {
  experiment = NULL;
  model = NULL;
}

SimulationBuilder::SimulationBuilder(std::string fn, bool verbose) {
  if (verbose){
    std::cout << "===== SETUP =====\n";
    std::cout << "Read simulation settings from " << fn << std::endl;
  }
  experiment = NULL;
  model = NULL;
  ParseFile(fn);
  experiment->model = model;
  if (verbose) {
    std::cout << "--- EXPERIMENT - " << experiment->type << " ---\n";
    experiment->PrintInfo();
    std::cout  << "--- MODEL ---\n";
    model->PrintInfo();
    std::cout << "===== INITIALIZATION =====\n";
  }
  std::map<std::string, std::string> pars = experiment->GetPars();
  for (auto item : pars){
    std::cout << item.first << " = " << item.second << std::endl;
  }

    model->gen.seed(experiment->init_seed);
  experiment->Initialize();
  model->gen.seed(experiment->seed);
  if (experiment->save_xml) {
    std::string xml_out = experiment->xml_out.empty() ? experiment->data_path + "/" + experiment->simulation_name + ".xml"
                                    : experiment->xml_out;
    this->ToXML(xml_out);
    experiment->out_files.insert(xml_out);
  }
}


void SimulationBuilder::ToXML(std::string fn) {
  tinyxml2::XMLDocument doc;
  tinyxml2::XMLNode *pRoot = doc.NewElement("Simulation");
  tinyxml2::XMLElement *base_el;

  time_t rawtime;
  char timestamp [80];
  time (&rawtime);
  struct tm * timeinfo = localtime (&rawtime);
  strftime (timestamp,80,"%Y-%m-%d %H:%M:%S",timeinfo);
  (pRoot->ToElement())->SetAttribute("run",timestamp);

  // Write model parameters
  base_el = doc.NewElement("Model");
  for (auto parset : model->GetPars()) {
    tinyxml2::XMLElement *el = doc.NewElement(parset.first.c_str());
    if (parset.second.size() > 0) { el->SetText(parset.second.c_str()); }
    base_el->InsertEndChild(el);
  }
  pRoot->InsertEndChild(base_el);

  // Write experiment parameters
  base_el = doc.NewElement("Experiment");
  base_el->SetAttribute("type", experiment->type.c_str());
  for (auto parset : experiment->GetPars()) {
    tinyxml2::XMLElement *el = doc.NewElement(parset.first.c_str());
    if (parset.second.size() > 0) { el->SetText(parset.second.c_str()); }
    base_el->InsertEndChild(el);
  }
  pRoot->InsertEndChild(base_el);
  doc.InsertFirstChild(pRoot);
  XMLError eResult = doc.SaveFile(fn.c_str());
}

std::map<std::string, std::string> SimulationBuilder::ParsePars(const tinyxml2::XMLNode *node) {
  std::map<std::string, std::string> pars;
  const tinyxml2::XMLElement *el = node->FirstChildElement();
  while (el) {
    if (el->FirstChild())
      pars[el->Name()] = el->GetText();
    else
      pars[el->Name()] = "";
    el = el->NextSiblingElement();
  }
  return pars;
}

void SimulationBuilder::BuildModel(const tinyxml2::XMLNode *node) {
  model = new GrowthModel(this->ParsePars(node));
}

void SimulationBuilder::BuildExperiment(const tinyxml2::XMLNode *node) {
  if ((node->ToElement())->Attribute("type", "IteratedGrowth")) {
    experiment = new ExperimentIteratedGrowth(this->ParsePars(node));
  }
  else if ((node->ToElement())->Attribute("type", "ConstantGrowth")) {
    experiment = new ExperimentConstantGrowth(this->ParsePars(node));
  }
}

void SimulationBuilder::ParseFile(std::string fn) {
  tinyxml2::XMLDocument xmlDoc;
  xmlDoc.LoadFile(fn.c_str());
  tinyxml2::XMLNode *root = xmlDoc.LastChild();
  const tinyxml2::XMLNode *node = root->FirstChild();
  std::string node_name;
  while (node) {
    node_name = std::string(node->Value());
    if (node_name.compare("Model") == 0)
      BuildModel(node);
    else if (node_name.compare("Experiment") == 0)
      BuildExperiment(node);
    node = node->NextSibling();
  }
}
