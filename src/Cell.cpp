#include "Cell.h"
#include <iostream>
#include <sstream>

Cell::Cell(){
  division_rate = 0;
  death_rate = 0;
  barcode = 0;
  dead = false;
  uid = 0;
  div_time = -1;
  death_time = -1;
  parent_uid = -1;
  birth_time = 0;
}

Cell::Cell(double _division_rate, double _death_rate, int _barcode, int _uid){
  division_rate = _division_rate;
  death_rate = _death_rate;
  barcode = _barcode;
  dead = false;
  uid = _uid;
  div_time = -1;
  death_time = -1;
}

void Cell::Print(){
  std::cout << "division rate = " << division_rate << std::endl;
  std::cout << "death rate = " << death_rate << std::endl;
  std::cout << "barcode = " << barcode << std::endl;
  std::cout << "dead = " << dead << std::endl;
  std::cout << "uid = " << uid << std::endl;
}



void Cell::Mutate(base_generator_type& gen){}

double Cell::GetDivTime(base_generator_type& gen){
  std::uniform_real_distribution<double> unif_distr(0.0,1.0);
  if (division_rate > 0)
    return -log(unif_distr(gen)) / division_rate;
  else
    return 0;
}

double Cell::GetDeathTime(base_generator_type& gen){
  std::uniform_real_distribution<double> unif_distr(0.0,1.0);
  if (death_rate > 0)
    return -log(unif_distr(gen)) / death_rate;
  else
    return 0;
}
