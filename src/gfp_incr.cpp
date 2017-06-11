#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/ref.hpp>
#include <models/fp_task.hpp>
#include <models/task_parser.hpp>
#include <analysis/hyperplane.hpp>
#include <analysis/fp_response_time.hpp>
#include <analysis/global_fp.hpp>
#include <analysis/global_fpp.hpp>
#include <analysis/global_edf.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/foreach.hpp>

#include <iomanip>
#include <map>

using namespace std;
using namespace Scan;

void print(vector<Task> & tasks) {
  for ( auto & x : tasks) {
    cout << x.get_wcet() << ", " << x.get_dline() << ", " << x.get_period() << endl;
  }
  cout << endl;
}

static int seed = 12345;
static int ul = 0;

Task gen_a_task(double Umin=0.05, double Umax=0.45, const int per = 0) {


  int Tmin = 1, Tmax = 200;

  //vector<int> ps = {10, 20, 50, 100};
  //vector<double> us = {0.05, 0.25, 0.45};
  //if ( ul ++ % 2 == 0) {Umin = 0.45; Umax = 0.5;}
  int p = 0, c = 0, d = 0;
  while ( true) {
    srand (seed ++);
    //p = ps[rand() % ps.size()]; //(Tmax - Tmin) + Tmin;
    p = rand() % (Tmax - Tmin) + Tmin;
    //if (per > 0 ) p = per;

    double dc = Umin;
    if (Umin != Umax) {
      default_random_engine generator(seed++);
      exponential_distribution<double> disC(0.1);
      //uniform_real_distribution<double> disC(Umin, Umax);
      //srand (seed ++);
      //dc = us[rand() % us.size()]; //disC(generator);
      dc = disC(generator);
      if (dc <= Umin) dc = Umin; //continue; 
      if (dc > Umax ) dc = Umax; //continue; 
    }
    c = dc * p;

    //default_random_engine generatorD(seed++);
    //uniform_real_distribution<double> disD(dc, 1);
    //double dd = disD(generatorD);
    //d =  dd * p; 
    d = p;
    if ( c > 0 and c < p and c <= d) break;
  }

  Task t(c, d, p);

  return t;

}

vector<Task> gen_a_taskset(const int m) {

  vector<Task> tasks;
  for ( int i = 0; i < m; i++)
    tasks.push_back(gen_a_task());
      
  return tasks;

}

double totU(vector<Task> & tasks) {
  double u = 0;
  for ( auto &x: tasks)
    u += x.get_wcet() * 1.0 / x.get_period();

  return u;

}


void gfp_run_tests_incr(const int m, const double xxx) {

  double fac = 0.1;

  int total = 1000000*1000;
  int num = 0;

  map<int, int> tot;
  map<int, int> rta;
  map<int, int> bl;

  for ( int i = 10*fac; i <= m * 10; i += 10*fac) {
    tot[i] = 0;
    rta[i] = 0;
    bl[i] = 0;
  }

  vector<Task> T = gen_a_taskset(m);

  while ( num ++ < total) {

    Task t = gen_a_task();
    T.push_back(t);

    // DM
    sort(T.begin(), T.end(), 
        [](const Task& a, const Task& b) {return a.get_dline() < b.get_dline();}
        );
      

    double util = totU(T);
    

    int step = (int) (util / fac); // * 10;
    if ( step >= m* 10) {
      T.clear();
      T = gen_a_taskset(m);
      continue;
    }

    bool brta = RTA_LC(T, m, true);
   // cout << brta << ", --- " << endl;
   // print(T);
    bool bbl = true; //baselinep(T, m);
    
    //cout << bbl << ", ++++ " << endl;
    //if ( brta and not bbl) print(T);
    if ( not brta) {
      bbl = baseline(T, m);
    }


    //{
    if ( tot[step] < 100000) {
      if (brta) rta[step] ++;
      if (bbl) bl[step] ++;
      tot[step] ++;
    }

    if ( not (brta or bbl)) {
      T.clear();
      T = gen_a_taskset(m);
    }

    if ( num % 10 == 0) {
      string fres = to_string(m) + ".txt-8"; //final-day-0.1-1-2000-0-1-constrD";
      //string fres = to_string(m) + ".txt-final-day-0.1-1-2000-0-1-constrD";
      ofstream of_res;
      of_res.open(fres);
      //of_res.open(fres, ios::app);


      for ( int i = 10*fac; i <= 10 * m; i += 10*fac)  {
        of_res << i << "  " ;
        of_res << rta[i] << "  " ;
        of_res << bl[i] << "  " ;
        of_res << tot[i] <<  endl ;
      }
      of_res << num <<  endl ;

      of_res.close();
    }


  }
  string fres = "final-" + to_string(m) + ".txt";
  ofstream of_res;
  of_res.open(fres, ios::app);


  for ( int i = 10*fac; i <= 10 * m; i += 10*fac)  {
    of_res << i << "  " ;
    of_res << rta[i] << "  " ;
    of_res << bl[i] << "  " ;
    of_res << tot[i] << endl ;
  }
  of_res << endl ;
  of_res.close();


}

int main(int argc, char *argv[])
{

  //Task t1 (10,20,20);
  //Task t2 (15, 30, 30);
  //Task t3 (25, 50, 50);
  //vector<Task> tasks = {t1, t2, t3};
  //bool rta = RTA_LC(tasks, 2, true);
  //bool bbl = baseline(tasks, 2);
  //cout << rta << endl << bbl << endl;

  try {
    gfp_run_tests_incr(8, 0.1);
  } catch (string s) {cout << s << endl;}
  ////gfp_run_tests_incr( 2, 0.1);
  //gfp_run_tests_incr( 10, 0.1);

  return 0;

}
