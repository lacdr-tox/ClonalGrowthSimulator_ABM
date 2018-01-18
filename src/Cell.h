#ifndef ABM_CELL_H
#define ABM_CELL_H

#include <random>
#include <string>
typedef	std::mt19937 base_generator_type;

class Cell {
 public:
  Cell();
  Cell(double _division_rate, double _death_rate, int _barcode, int _uid);
  void Print();

  void Mutate(base_generator_type & gen);

  double GetDivTime(base_generator_type &gen);

  double GetDeathTime(base_generator_type &gen);

  bool dead;
  double division_rate, death_rate;
  int barcode, uid;
  double div_time, death_time, birth_time;
  int parent_uid;

};


#endif //ABM_CELL_H
