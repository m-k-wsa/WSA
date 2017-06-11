#include "task.hpp"
#include "wa.hpp"
#include "exp.hpp"
#include <time.h>

#include <iostream>
using namespace std;

const int bigN=1000000000;

int main(int argc, char* argv[])
{
  if(argc!=3)
  {
    std::cout << "Please provide m and K as the input arguments.\n";
    exit(1);
  }
  int e=1000;
  Task t1 (e+15.83, e+15.83, 20000, 20000);
  t1.cs=0.41;
  Task t2(2309.5, 2309.5, 20000, 20000);
  t2.cs=105.42;
  Task t3(1148.64, 1148.64, 25000, 25000);
  t3.cs=41.41;
  Task t4(2016.33*30/25, 2016.33*30/25, 30000, 30000);
  t4.cs=29.75;
  Task t5(287.5, 287.5, 50000, 50000);
  t5.cs=247.06;
  Task t6(42.56*6/5, 42.56*6/5, 60000, 60000);
  t6.cs=7.77;
  Task t7(2318.42, 2318.42, 100000, 100000);
  t7.cs=237.98;
  Task t8(5296.84, 5296.84, 100000, 100000);
  t8.cs=137.13;
  Task t9(325.64, 325.64, 200000, 200000);
  t9.cs=14.97;
  Task t10(3285.24, 3285.24, 200000, 200000);
  t10.cs=240.13;
  Task t11(208.67, 208.67, 500000, 500000);
  t11.cs=72.22;
  Task t12(539.5, 539.5, 500000, 500000);
  t12.cs=35.01;
  Task t13(47616.35, 47616.35, 1000000, 1000000);
  t13.cs=0;
  Task t14(799005.74, 799005.74, 2000000, 2000000);
  t14.cs=318.37*1000;
  Task t15(e*e+406.02, e*e+406.02, 10000000, 10000000);
  t15.cs=0;

  vector<Task> tasks={t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15};

  double util=0;
  for(int i=0; i<tasks.size(); i++)
  {
    util+=tasks[i].wcet/tasks[i].period;
    cout << "Total util: " << util << "\n";
    vector<Task> hps;
    for(int j=0; j<i; j++)
      hps.push_back(tasks[j]);
    compute_wcrt(tasks[i], hps);
    compute_bp(tasks[i], hps);
    compute_bcrt(tasks[i], hps);
    tasks[i].print();
    cout << "Total util: " << util << "\n";

  }
  vector<double> CS={318.37*1000, 12.27, 12.27, 7.26, 7.26, 7.26};
  int nRes=1;
  int i=14;
  vector<Task> hps, lps;
  for(int j=1; j<=15; j++)
  {
    if(j<i)
    {
      hps.push_back(tasks[j-1]);
      hps[j-1].print();
    }
    if(j>i)lps.push_back(tasks[j-1]);
  }
  // this is a shortcut 
  hps[0].CS[0]=CS;
  lps[0].CS[0]=CS;
  int m=atoi(argv[1]), K=atoi(argv[2]);
  tasks[i-1].print();
  statWA swa = wanalysis(tasks[i-1], hps, m, K, lps, nRes);
  cout << "hps size: " << hps.size() << ", (K, m) " << K << ", " << m << ", bounded = " << swa.bounded << ", time = " << swa.time << ", R=" << tasks[i-1].wcrt<< ", bp=" << tasks[i-1].bp << endl;

  return 0;
}
