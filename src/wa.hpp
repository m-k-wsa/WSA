#ifndef __WA_HPP__
#define __WA_HPP__

#include "task.hpp"
#include <vector>

struct statWA {
  double time;
  bool bounded;
};



void hello_cplex();

//int wanalysis(const Task& ti, const int K, const std::vector<Task>& hps);
statWA wanalysis(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

#endif
