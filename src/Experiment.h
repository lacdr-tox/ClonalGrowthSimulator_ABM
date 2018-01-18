#ifndef ABM_EXPERIMENT_H
#define ABM_EXPERIMENT_H

#include "GrowthModel.h"
#include <boost/random/mersenne_twister.hpp>
#include <cinttypes>
#include <set>

class Experiment {

 public:
  Experiment();

  Experiment(std::map<std::string, std::string> pars);

  void SetPars(std::map<std::string, std::string> pars);
  std::map<std::string, std::string> GetBasePars();

  virtual ~Experiment();

  virtual void Run();

  std::vector<int> GetBarcodesUniform(int nbarcodes);
  std::vector<int> GetBarcodesUniform(int nbarcodes, int nof_cells);
  std::vector<int> GetBarcodesUniform(int nbarcodes, double ratio, int nof_cells);
  std::vector<int> GetBarcodesFromFile(std::string init_clones_file);
  std::vector<int> GetBarcodesFromFile(std::string init_clones_file, int nof_cells);
  std::vector<int> GetBarcodesFromFile(std::string init_clones_file, double ratio, int nof_cells);
  std::vector<int> GetNoise();

  void GrowToSize(int ncells_max);
  void GrowForTime(double tmax);
  void SelectCells(int ncells);
  void SelectBarcodedCells();
  //!
  /*!
   * Initialize clones based on the values of class parameters \p init_from_file, \p init_clones_file,
   * \p ncells_init, and nclones.
   *
   */
  void Initialize();

  //!
  /*!
   * Print some information about the experiment setup
   *
   */
  void PrintInfo();

  //!
  /*!
   * Get map with all experiment parameters.
   *
   */
  virtual std::map<std::string, std::string> GetPars() { };

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   */
  void SaveState();

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   *
   * \param time time step to be put in the first column
   */
  void SaveState(double time);

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   *
   * \param time time step to be put in the first column
   * \param postfix of the output file
   */
  void SaveState(double time, std::string postfix);
  void SaveDivRates(double time);
  void SaveDivRatesBC(double time);
  void SaveDivRates(double time, std::string postfix);
  void SaveDivRatesBC(double time, std::string postfix);
  void SavePopulationStats(double time);
  void SaveCellInfo(std::string postfix);
  //!
  /*!
   * Gzip all output files
   *
   */
  void GzipOutput();

  //!
  /*!
   * Set up folder for simulation results
   *
   */
  void SetupSimulationFolder();

  //! Seed for simulation algorithm (growth and passage)
  unsigned seed;
  //! Seed for initialization
  unsigned init_seed;
  unsigned init_seed_distr;
  bool init_from_pool_file;
  std::string pool_file;
  //! Filename of file with initial clone distribtion
  std::string init_clones_file;
  //! Initialize clones using a stored distribution
  bool init_from_file,init_large_pool;
  //! Number of cells to initialize
  unsigned int ncells_init;
  //! Growth model that contains all possible transitions
  GrowthModel *model;
  //! Experiment type
  std::string type;
  //! Simulation identifier
  std::string simulation_name;
  //! Path to output data
  std::string data_path;
  //! Gzip generated data
  bool gzip;
  //! List of output files
  std::set <std::string> out_files;
  std::string xml_out;
  bool save_xml;

 protected:
  //!
  /*!
   * Read initial clone distribution from file
   *
   * \param init_clones_file name of file with input data
   * \return clone fraction (cnt/n_reads) per clone
   *
   */
  std::vector<int> ReadInitialClones(std::string init_clones_file);
  std::vector<double> ReadInitialCloneFreq(std::string init_clones_file);
  std::vector< std::pair<int,double>> ReadInitPool(std::string fn);

  //! Store results in separate folder
  bool use_simdir;
  //! Store simulation results
  bool write_output;
  //! Store population size
  bool write_pop_size;
  //! Sum of all clone frequencies
  std::uint64_t population_size;
  //! Current simulation step
  unsigned int current_step;
  int nclones;
  bool track_divtime, track_divtime_bc,track_cells;
  double infection_ratio;
  int ncells_infection;
  double post_infection_growth_time,post_selection_growth_time;
  double current_time;
  bool add_noise;
  int min_size_noise,max_size_noise;
  int nof_noise;

 private:
  bool is_number(const std::string& s);
  bool monoclonal_full,monoclonal_simple;

};


#endif //ABM_EXPERIMENT_H
