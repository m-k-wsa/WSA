#include "task.hpp"
#include "wa.hpp"
#include "exp.hpp"
#include <time.h>

#include <iostream>
using namespace std;

int main() 
{

  double jitter=20*1000; //20
  Task t1(jitter,40*1000, 40*1000, 100*1000, 100*1000, 100*1000); 
  Task t2(0, 17*1000, 17*1000, 40*1000, 39999, 40001); 
  t1.bcrt=t1.wcet;
  t1.wcrt=t1.wcet+2000;
  compute_wcrt(t2, {t1});
  cout << "t2 wcrt: " << t2.wcrt << "\n";
  compute_bp(t2, {t1});
  cout << "t2 bp: " << t2.bp << "\n";
  t2.bcrt=t2.wcet;

  int m, K;
  statWA res = 
  //wanalysis(t2, {t1}, 40, 100);
  //wanalysis(t2, {t1}, 6, 15);
  //wanalysis(t2, {t1}, 3, 10);
  wanalysis(t2, {t1}, 2, 5);
  //wanalysis(t2, {t1}, 2, 4);
  //wanalysis(t2, {t1}, 2, 3);
  //wanalysis(t2, {t1}, 1, 2);
  //wanalysis(t2, {t1}, 1, 3);
  std::cout << "bounded = " << res.bounded << "\n"; 

  return 0;
}
