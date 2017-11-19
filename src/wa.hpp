#ifndef __WA_HPP__
#define __WA_HPP__
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "task.hpp"
#include <vector>

struct statWA {
  double time;
  bool bounded;
};



void hello_cplex();

std::string  convert_to_string(const int Number);
//int wanalysis(const Task& ti, const int K, const std::vector<Task>& hps);
statWA wanalysis(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

statWA wanalysis_kill(const Task& ti, const std::vector<Task>& hps, const int m, const int K);

//void wsa(std::vector<Task> &tasks, m, K);

#endif
