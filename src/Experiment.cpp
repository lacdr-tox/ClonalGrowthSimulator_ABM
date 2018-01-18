#include "Experiment.h"

#include <iostream>
#include <fstream>      // std::ofstream
#include "Experiment.h"
#include <chrono>
#include <algorithm>    // std::min_element, std::max_element
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/filesystem.hpp>
#include <limits>       // std::numeric_limits
#include <utility>      // std::pair, std::make_pair


Experiment::Experiment() {
  seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
  init_seed = 0;
  init_seed_distr = 0;
  current_step = 0;
  gzip = false;
  use_simdir = false;
  write_output = false;
  write_pop_size = false;
  data_path = ".";
  type = "";
  simulation_name = "";
  init_clones_file = "";
  init_from_file = false;
  ncells_init = 0;
  nclones = 1;
  track_divtime = false;
  track_divtime_bc = false;
  save_xml = false;
  current_time = 0;
  add_noise = false;
  min_size_noise = 1;
  max_size_noise = 1;
  nof_noise = 0;
  track_cells = false;
  init_large_pool = true;
  init_from_pool_file = false;
}

Experiment::Experiment(std::map<std::string, std::string> pars) : Experiment(){
  SetPars(pars);
}


void Experiment::SetPars(std::map<std::string, std::string> pars){
  for (auto par : pars) {
    if (par.first.compare("Name") == 0) {
      simulation_name = par.second;
      write_output = true;
    }
    else if (par.first.compare("WritePopSize") == 0)
      write_pop_size = true;
    else if (par.first.compare("OutPath") == 0)
      data_path = par.second;
    else if (par.first.compare("InitialPopulationSize") == 0)
      ncells_init = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitFile") == 0) {
      init_clones_file = par.second;
      init_from_file = true;
    }
    else if (par.first.compare("InitUniform") == 0) {
      nclones = (unsigned int) std::stoi(par.second);
      init_from_file = false;
    }
    else if (par.first.compare("InitFromPoolFile") == 0) {
      init_from_pool_file = true;
      init_large_pool = false;
      pool_file = par.second;
    }
    else if (par.first.compare("InitSimple") == 0)
      init_large_pool = false;
    else if (par.first.compare("Seed") == 0)
      seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitSeed") == 0)
      init_seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitSeedDistr") == 0)
      init_seed_distr = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("gzip") == 0)
      gzip = true;
    else if (par.first.compare("UseSimDir") == 0)
      use_simdir = true;
    else if (par.first.compare("TrackDivTimeBC") == 0)
      track_divtime_bc = true;
    else if (par.first.compare("TrackDivTime") == 0)
      track_divtime = true;
    else if (par.first.compare("TrackCells") == 0)
      track_cells = true;
    else if (par.first.compare("SaveXML") == 0){
      save_xml = true;
      xml_out = par.second;
    }
    else if (par.first.compare("Monoclonal") == 0)
      monoclonal_simple = true;
    else if (par.first.compare("MonoclonalSimple") == 0)
      monoclonal_simple = true;
    else if (par.first.compare("MonoclonalFull") == 0)
      monoclonal_full = true;
    else if (par.first.compare("InfectionRatio") == 0)
      infection_ratio = std::stod(par.second);
    else if (par.first.compare("InfectionPopulationCells") == 0)
      ncells_infection = std::stoi(par.second);
    else if (par.first.compare("PostInfectionGrowthTime") == 0)
      post_infection_growth_time = std::stod(par.second);
    else if (par.first.compare("PostSelectionGrowthTime") == 0)
      post_selection_growth_time = std::stod(par.second);
    else if (par.first.compare("AddNoise") == 0){
      add_noise = true;
      if (par.second.size() > 0)
        nof_noise = std::stoi(par.second);
    }
    else if (par.first.compare("NoiseMinSize") == 0)
      min_size_noise = std::stoi(par.second);
    else if (par.first.compare("NoiseMaxSize") == 0)
      max_size_noise = std::stoi(par.second);
  }
  if (init_seed == 0) { init_seed = seed; }
  if (init_seed_distr == 0){ init_seed_distr = init_seed;}
  if (use_simdir){ SetupSimulationFolder(); }
}

std::map<std::string, std::string> Experiment::GetBasePars(){
  std::map<std::string, std::string> pars = {{"OutPath", data_path},
                                             {"InitialPopulationSize", std::to_string(ncells_init)},
                                             {"Seed", std::to_string(seed)},
                                             {"InitSeed", std::to_string(init_seed)}};
  if (monoclonal_simple) { pars["MonoclonalSimple"]= "";}
  if (monoclonal_full) {
    pars["MonoclonalFull"]= "";
    pars["InfectionRatio"] = std::to_string(infection_ratio);
    pars["InfectionPopulationCells"] = std::to_string(ncells_infection);
    pars["PostInfectionGrowthTime"] = std::to_string(post_infection_growth_time);
    pars["PostSelectionGrowthTime"] = std::to_string(post_selection_growth_time);
  }
  if (add_noise){
    pars["AddNoise"] = std::to_string(nof_noise);
    pars["NoiseMinSize"] = std::to_string(min_size_noise);
    pars["NoiseMaxSize"] = std::to_string(max_size_noise);
  }
  if (save_xml) { pars["SaveXML"] = xml_out; }
  if (gzip) { pars["gzip"] = ""; }
  if (!init_from_file){ pars["InitUniform"] = std::to_string(nclones);}
  else{ pars["InitFile"] = init_clones_file;}
  if (write_output){ pars["Name"] = simulation_name;}
  if (write_pop_size) { pars["WritePopSize"] = "";}
  if (use_simdir) { pars["UseSimDir"] = ""; }
  if (track_divtime_bc){ pars["TrackDivTimeBC"] = "";}
  if (track_divtime){ pars["TrackDivTime"] = "";}
  if (init_from_pool_file){ pars["InitFromPoolFile"] = pool_file;}
  if (!init_large_pool){ pars["InitSimple;"] = "";}
  return pars;
};


void Experiment::PrintInfo(){
  std::map<std::string, std::string> pars = this->GetPars();
  for (auto par : pars){
    if (par.second.size() > 0)
      std::cout << par.first << " = " << par.second << std::endl;
    else
      std::cout << par.first << " = True" << std::endl;
  }
};

Experiment::~Experiment() { }

void Experiment::SaveCellInfo(std::string postfix){
  std::ofstream of;
  std::string fn = data_path + "/" + simulation_name + "_cells_" + postfix + ".txt";
  of.open(fn, std::ios::trunc);
  out_files.insert(fn);
  of << "uid\tbarcode\tdivision_rate\tdeath_rate";
  for (auto item : model->cells){
    of << "\n" << item.second.uid << "\t" << item.second.barcode << "\t";
    of << item.second.division_rate << "\t" << item.second.death_rate;
  }
  of.close();
}

void Experiment::SaveState() { SaveState(current_step); }

void Experiment::SaveState(double time) {
  SaveState(time, "clones");
}

void Experiment::SaveState(double time, std::string postfix) {
  std::ofstream of;
  std::string fn = data_path + "/" + simulation_name + "_" + postfix + ".txt";
  if (out_files.find(fn) == out_files.end()){
    of.open(fn, std::ios::trunc);
    out_files.insert(fn);
  }
  else {of.open(fn, std::ios::app);}
  std::string s = "\n" + std::to_string(time);
  for (auto item : model->barcodes){s += "\t" + std::to_string(item.second); }
  of << s;
  of.close();
}

void Experiment::SaveDivRatesBC(double time) {
  SaveDivRatesBC(time,"divratesBC");
}

void Experiment::SaveDivRatesBC(double time, std::string postfix) {
  std::ofstream of;
  std::string fn = data_path + "/" + simulation_name + "_" + postfix + ".txt";
  if (out_files.find(fn) == out_files.end()){
    of.open(fn, std::ios::trunc);
    out_files.insert(fn);
    of << time;
  }
  else {
    of.open(fn, std::ios::app);
    of << "\n" << time;
  }
  std::map <int, double> mean_div_rate;
  for (auto item : model->cells){
    mean_div_rate[item.second.barcode] += item.second.division_rate/model->barcodes[item.second.barcode];
  }
  for (auto item : mean_div_rate){ of << "\t" << item.second; }
  of.close();
}

void Experiment::SaveDivRates(double time) {
  SaveDivRates(time,"divrates");
}

void Experiment::SaveDivRates(double time,std::string postfix) {
  std::ofstream of;
  std::string fn = data_path + "/" + simulation_name + "_" + postfix + ".txt";
  if (out_files.find(fn) == out_files.end()){
    of.open(fn, std::ios::trunc);
    out_files.insert(fn);
    of << time;
  }
  else {
    of.open(fn, std::ios::app);
    of << "\n" << time;
  }
  for (auto item : model->cells){
    of << "\t" << item.second.division_rate;
  }
  of.close();
}

void Experiment::GzipOutput() {
  for (auto f : out_files){
    std::cout << "gzip " << f << std::endl;
    std::ifstream inStream(f, std::ios_base::in);
    std::ofstream outStream(f+".gz", std::ios_base::out);
    boost::iostreams::filtering_streambuf< boost::iostreams::input> in;
    in.push( boost::iostreams::gzip_compressor());
    in.push( inStream );
    boost::iostreams::copy(in, outStream);
    remove(f.c_str());
  }
}



void Experiment::SavePopulationStats(double time){
  // Create list of frequencies for the population
  double fleft,f50,gini,fleft_real;
  std::vector <uint64_t> freq(model->barcodes.size(),0);
  for (auto bc : model->barcodes){freq.push_back(bc.second);}
  std::cout << freq.size() << " " << current_step << " " << model->nclones0 << std::endl;

  if (freq.size() > 0) {

    // compute the fraction of clones that is left, compared to the inital state
    unsigned int nleft = 0;
    for (auto cnt : freq) { if (cnt > 0) { nleft++; }}
    if (current_step == 0) { model->nclones0 = nleft; }
    fleft = nleft / ((double) model->nclones0);
    fleft_real = nleft / ((double) model->real_nclones0);
    // compute fraction of clones that makes up the first 50 percent of the population
    std::vector<uint64_t> freq_nz;
    for (auto cnt : freq) { if (cnt > 0) { freq_nz.push_back(cnt); }}
    std::sort(freq_nz.begin(), freq_nz.end(), std::greater<int>());
    std::vector<uint64_t> sfreq(nleft, 0);
    sfreq[0] = freq_nz[0];
    for (unsigned int i = 1; i < nleft; i++) { sfreq[i] = sfreq[i - 1] + freq_nz[i]; }
    f50 = 100 * sfreq[(int) (.5 * nleft)] / ((double) model->nof_cells);

    // compute gini coefficient
    std::vector <double> frac;
    for (auto cnt : freq_nz){frac.push_back(cnt/((double)model->nof_cells));}
    double sdiff = 0;
    for (auto f : frac){
      for (auto ff : frac){ sdiff += fabs(f-ff); }
    }
    gini = sdiff/(2.*nleft);
  }
  else{
    fleft = 1;
    fleft_real = 1;
    f50 = 0;
    gini = 0;
  }

  // compute division rate mean
  double mu_divrate = 0;
  for (auto item : model->cells){ mu_divrate += item.second.division_rate/model->nof_cells;}

  // write results to file
  std::ofstream of;
  std::string fn = data_path + "/" + simulation_name + "_popstats.txt";
  if (out_files.find(fn) == out_files.end()){
    of.open(fn, std::ios::trunc);
    of << "#time\tpopsize\tfleft\tf50\tgini\tmu_divrate\tfleft_real";
    out_files.insert(fn);
  }
  else {of.open(fn, std::ios::app);}

  of << "\n" << time << "\t" << model->nof_cells << "\t" << fleft << "\t" << f50 << "\t" << gini;
  of << "\t" << mu_divrate << "\t" << fleft_real;
  of.close();
}


void Experiment::SetupSimulationFolder(){
  std::string data_path_new = data_path + "/" + simulation_name + "/";
  boost::filesystem::path p (data_path_new);
  // clean up old simulation files
  boost::filesystem::create_directory(p);
  data_path = data_path_new;
}

void Experiment::Run() {
  std::cout << "Run undefined experiment\n";
}

void Experiment::Initialize(){
  int ncells;
  current_time = 0;
  if (monoclonal_full){
    std::cout << "Perform full monoclonal initialization\n";
    model->AddCells(1);
    SavePopulationStats(current_time);
    this->GrowToSize(ncells_infection);
    std::vector<int> barcode_distribution;
    if (init_from_file)
      barcode_distribution= GetBarcodesFromFile(init_clones_file,infection_ratio);
    else
      barcode_distribution = GetBarcodesUniform(nclones,infection_ratio);
    model->AssignBarcodes(barcode_distribution);

    // Note that we first select and then grow.
    // Since we assume that cells do not affect the growth of other cells
    // the order of selection and growth cannot affect the outcome.
    // However, by switching the order, we do not have to simulate cells
    // that will be discarded later.
    this->SelectBarcodedCells();
    this->GrowForTime(post_infection_growth_time+post_selection_growth_time);
    this->SelectCells(ncells_init);
    SaveDivRates(current_time,"divrates_init");
  }
  else if (init_large_pool) {
    std::vector<int> init_clones;
    if (init_from_file) {
      init_clones = ReadInitialClones(init_clones_file);
    }
    else{
      int cells_per_clone = (int)(ceil((double)ncells_init/nclones));
      for (int i = 0; i < nclones; i++){init_clones.push_back(cells_per_clone);}
    }
    int uid = 0;
    std::mt19937 gen(init_seed);
    // create init population with all barcodes
    std::normal_distribution<double> distr(1, model->division_rate_sd);
    model->real_nclones0 = init_clones.size();
    for (unsigned int i = 0; i < init_clones.size(); i++) {
      for (int j = 0; j < init_clones[i]; j++) {
        double r = distr(gen) * model->division_rate;
        model->AddInitCell(uid, i, r);
        uid++;
      }
    }
    if (track_cells) { SaveCellInfo("distribution"); }
    model->max_uid = uid;
    model->gen.seed(init_seed_distr);
    std::vector<int> noise;
    if (add_noise) {
      noise = GetNoise();
      ncells_init = ncells_init- std::accumulate(noise.begin(), noise.end(), 0);
    }
    std::cout << ncells_init << std::endl;
    this->SelectCells(ncells_init);
  }
  else if (init_from_pool_file){
    std::cout << "Init from pool file " << pool_file << std::endl;
    std::mt19937 gen(init_seed);
    std::vector< std::pair<int,double> > pool = ReadInitPool(pool_file);
    std::vector<int> idx;
    for (int i = 0; i < (int)pool.size(); i++){idx.push_back(i);}
    std::shuffle(idx.begin(), idx.end(), gen);
    int uid = 0;
    for (unsigned int i = 0; i < ncells_init; i++){
      model->AddInitCell(uid, pool[i].first, pool[i].second);
      uid++;
    }
    model->max_uid = uid;
    model->nof_cells = ncells_init;
  }
  else{
    std::cout << "add " << ncells_init << " cells\n";
    model->AddCells(ncells_init);
    std::cout << model->nof_cells << std::endl;
    std::vector<int> barcode_distribution;
    if (init_from_file)
      barcode_distribution= GetBarcodesFromFile(init_clones_file,ncells_init);
    else
      barcode_distribution = GetBarcodesUniform(ncells_init);
    model->AssignBarcodes(barcode_distribution);
  }
  if (track_cells){SaveCellInfo("init");}
}

void Experiment::SelectCells(int ncells){
  std::cout << "select " << ncells << " out of " << model->cells.size() << "cells\n";
  std::vector<int> idx;
  for (auto c : model->cells){ idx.push_back(c.first);}
  std::shuffle(idx.begin(), idx.end(), model->gen);
  std::vector<int> keep(idx.begin(), idx.begin()+ncells);
  model->UpdateQueues(keep);
}

void Experiment::SelectBarcodedCells(){
  std::vector<int> keep;
  for (auto c : model->cells){
    if (c.second.barcode >= 0){ keep.push_back(c.first);}
  }
  std::cout << "Keep " << keep.size() << " barcoded cells out of " << model->cells.size() << " cells\n";
  model->UpdateQueues(keep);
}

void Experiment::GrowToSize(int ncells_max) {
  std::cout << "Grow from " << model->nof_cells << " to " << ncells_max << " cells\n";
  int p_next = 0;
  double s_next = .1;
  while (model->nof_cells < ncells_max){
    current_time = model->Step();
    if ((100 * model->nof_cells) / ncells_max > p_next) {
      std::cout << p_next << "% - " << model->nof_cells << " cells\n";
      p_next += 2;
    }
    if (current_time >= s_next){
      SavePopulationStats(current_time);
      s_next += .1;
    }
  }
}

void Experiment::GrowForTime(double tmax) {
  int p_next = 0;
  double t0 = current_time;
  std::cout << "Grow from " << t0 << " to " << t0+tmax << std::endl;
  while (current_time < (tmax+t0)){
    current_time = model->Step();
    if ((100*(current_time-t0))/tmax > p_next){
      std::cout << p_next << "% - " << model->nof_cells << " cells\n";
      p_next += 2;
      SavePopulationStats(current_time);
    }
  }
}

std::vector<int> Experiment::GetBarcodesFromFile(std::string init_clones_file){
  return GetBarcodesFromFile(init_clones_file, 1, model->nof_cells);
}

std::vector<int> Experiment::GetBarcodesFromFile(std::string init_clones_file, int nof_cells){
  return GetBarcodesFromFile(init_clones_file, 1, nof_cells);
}

std::vector<int> Experiment::GetBarcodesFromFile(std::string init_clones_file, double ratio, int nof_cells){
  std::mt19937 gen(init_seed_distr);
  std::cout << "Read barcodes from " << init_clones_file << std::endl;
  int ncells = (int)round(ratio*nof_cells);
  std::cout << nof_cells << " " << ratio << " " << ncells << std::endl;
  std::vector<double> init_clones = ReadInitialCloneFreq(init_clones_file);
  std::vector<int> barcode_distribution(init_clones.size(),0);
  boost::random::discrete_distribution<> dist(init_clones);
  for (int i = 0; i < ncells; i++)
    barcode_distribution[dist(gen)] += 1;
  return barcode_distribution;
}



std::vector<int> Experiment::GetBarcodesUniform(int nbarcodes) {
  return GetBarcodesUniform(nbarcodes, 1, model->nof_cells);
}

std::vector<int> Experiment::GetBarcodesUniform(int nbarcodes, int nof_cells) {
  return GetBarcodesUniform(nbarcodes, 1, nof_cells);
}

std::vector<int> Experiment::GetBarcodesUniform(int nbarcodes, double ratio, int nof_cells){
  std::mt19937 gen(init_seed_distr);
  std::cout << "Uniformly assign " << nbarcodes << " barcodes\n";
  int ncells = (int)round(ratio*nof_cells);
  std::cout << "initialize uniform\n";
  int nperbc = (int)floor((double)ncells)/((double)nbarcodes);
  std::vector<int> barcode_distribution (nbarcodes,nperbc);
  int bd_count = nbarcodes*nperbc;
  std::vector<int> idx;
  for (int i = 0; i < nbarcodes; i++) { idx.push_back(i); }
  std::shuffle(idx.begin(), idx.end(), gen);
  for (int i = 0; i < (ncells-bd_count); i++){ barcode_distribution[idx[i]] += 1; }
  return barcode_distribution;
}

std::vector<int> Experiment::GetNoise(){
  std::uniform_int_distribution<int> distr(min_size_noise,max_size_noise);
  std::vector<int> noise;
  for(int i = 0; i < nof_noise; i++){
    noise.push_back(distr(model->gen));
  }
  return noise;
};

std::vector<int> Experiment::ReadInitialClones(std::string fn) {
  std::ifstream init_file(fn);
  std::string Barcode;
  std::string CloneSize;
  std::vector<int> reads;
  long double nof_reads = 0;
  while (init_file >> Barcode >> CloneSize) {
    if (is_number(CloneSize)) {
      reads.push_back(atoi(CloneSize.c_str()));
      nof_reads += atoi(CloneSize.c_str());
    }
  }
  return reads;
}

std::vector<double> Experiment::ReadInitialCloneFreq(std::string fn) {
  std::ifstream init_file(fn);
  std::string Barcode;
  std::string CloneSize;
  std::vector<int> reads;
  std::vector<double> init_clones;
  long double nof_reads = 0;
  while (init_file >> Barcode >> CloneSize) {
    if (is_number(CloneSize)) {
      reads.push_back(atoi(CloneSize.c_str()));
      nof_reads += atoi(CloneSize.c_str());
    }
  }
  for (int cnt : reads) {
    init_clones.push_back(cnt / nof_reads);
  }
  return init_clones;
}


std::vector< std::pair<int,double> > Experiment::ReadInitPool(std::string fn) {
  std::ifstream file(fn, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(file);
  std::istream instream(&in);
  std::string line;
  std::string s_uid,s_barcode,s_divrate,s_deathrate;
  std::vector< std::pair<int,double> > pool;
  while (instream >> s_uid >> s_barcode >> s_divrate >> s_deathrate){
    if (!is_number(s_uid))
      continue;
    pool.push_back(std::make_pair(atoi(s_barcode.c_str()),atof(s_divrate.c_str())));
  }
  return pool;
}

bool Experiment::is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(),
                                    s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}
