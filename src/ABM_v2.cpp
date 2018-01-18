#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>  // for high_resolution_clock
#include <ctime>
#include <map>

typedef	 std::mt19937 base_generator_type;

struct Cell{
  double division_rate;
  int barcode;
};


int main() {
  base_generator_type gen (1);
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  int nbc = 1;
  int n0 = 3e5;
  int nmax = 4e6;
  double r = 1;

  // initialize
  std::multimap<double, Cell> div_queue;
  double divt;
  for (int i = 0; i < n0; i++){
    divt = -log(distribution(gen))/r;
    div_queue.insert(std::pair<double, Cell>(divt,{r,1}));
  };

  double t;
  double savefreq = 0.01;
  double savenext = 0;

  int p_next = 0;
  std::cout << "grow\n";
  std::ofstream of;
  std::string fn = "growth_v4.txt";
  of.open(fn, std::ios::trunc);
  of << "#time\tnumber of cells\tsimulation time";
  of.close();
  std::ofstream oft;
  oft.open("dt.txt",std::ios::trunc);
  oft << "#dt";
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "Growth from " << n0 << " to " << nmax << " cells\n";
  while (div_queue.size() < nmax){
    t = div_queue.begin()->first;
    oft << "\n" << t;
    div_queue.erase(div_queue.begin());
    div_queue.insert(std::pair<double, Cell>(t-log(distribution(gen))/r,{r,1}));
    div_queue.insert(std::pair<double, Cell>(t-log(distribution(gen))/r,{r,1}));
    if (t >= savenext){
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      of.open(fn, std::ios::app);
      of << "\n" << t << "\t" << div_queue.size() << "\t" << elapsed.count();
      of.close();
      savenext = savenext+savefreq > t ? savenext+savefreq : t+savefreq;
    }
    if (100*((double)div_queue.size())/nmax > p_next){
      std::cout << p_next << "%\n";
      p_next += 10;
    }
  }
  oft.close();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "simulation took " << elapsed.count() << " seconds" << std::endl;

  return 0;
}