
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "exp.hpp"
#include "util.hpp"

using namespace std;

static int seed = 123;

//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
//  std::stringstream ss(s);
//  std::string item;
//  while (std::getline(ss, item, delim)) {
//    elems.push_back(item);
//  }
//  return elems;
//}
//
//
//std::vector<std::string> split(const std::string &s, char delim) {
//  std::vector<std::string> elems;
//  split(s, delim, elems);
//  return elems;
//}

//bool com_dline(const Task& a, const Task& b) {return a.get_dline() < b.get_dline();}

void simul(const int n, const int nSets, const string& fname, const int m, const int K, const int lb, const int ub) {

  ifstream input(fname);
  int count = 0;
  int nsched = 0;

  while (count ++ < nSets) {
    if(count>=ub) break; 

    vector<Task> tasks;
    string line;

    bool valid=true;
    for ( int i = 0; i < n; i ++) {
      getline(input, line); // read a line from the input file
      stringstream linestream(line);
      string data;
      getline(linestream, data, '\n');

      vector<string> tokens = split(data, ' ');

      double u;
      double c, p;
      u = atof(tokens[0].c_str());
      c = atof(tokens[1].c_str());
      p = atof(tokens[2].c_str());

      if(c==0 or p==0) valid=false;
      Task t(c, c, p, p);
      tasks.push_back(t);
    }
    getline(input, line); // between two tasksets, there is a line of delimiter
    
    if(count<lb) continue;

    if(not valid) continue;

    //tasks_sorting(tasks);
    sort(tasks.begin(), tasks.end(), //com_dline);
        [](const Task& a, const Task& b) {return a.get_dline() < b.get_dline();}
    );
    std::cout << "\n----------------------\n"; 
    for ( int x = 1; x <= n; x++) 
      tasks[x-1].print();


    // schedulability analysis
    for ( int x = 1; x <= n; x++) {
      cout << "x: " << x << "\n ";
      vector<Task> hps;
      for ( int y = 1; y < x; y++) 
        hps.push_back( tasks[y-1]);
      compute_wcrt(tasks[x-1], hps);
      //cout << "wcrt: " << tasks[x-1].get_wcrt() << ", " << tasks[x-1].get_dline() << "\n ";
      compute_bp(tasks[x-1], hps);
      //cout << "bp: " << tasks[x-1].get_bp() << "\n ";
      //tasks[x-1].print();
      compute_bcrt(tasks[x-1], hps);
      //cout << "bcrt: " << tasks[x-1].get_bcrt() << "\n ";
      if (x < n) 
      {
        //if(tasks[x-1].wcrt<=tasks[x-1].dline) continue;
        //else break;
        continue;
      }
      if(tasks[x-1].wcrt<=tasks[x-1].dline) continue; //{ std::cout << "sched: " << count << std::endl; continue;}
      if(tasks[x-1].bcrt>tasks[x-1].dline) continue;

      nsched ++;
      cout << "nsched, count: " << nsched << ", " << count << std::endl;
      //if( nsched > 100) return;

      for ( int y = 1; y < x; y++) 
        std::cout << "wcet = " << tasks[y-1].wcet << ", bcrt = " << tasks[y-1].get_bcrt() << ", wcrt = " << tasks[y-1].get_wcrt() << ", dline= " << tasks[y-1].get_dline() << ", bp=" << tasks[y-1].get_bp() << std::endl;
      std::cout << "wcet = " << tasks[x-1].wcet << ", bcrt = " << tasks[x-1].get_bcrt() << ", wcrt = " << tasks[x-1].get_wcrt() << ", dline= " << tasks[x-1].get_dline() << ", bp=" << tasks[x-1].get_bp() << std::endl;
      cout << "----------*--------------" << endl;

      statWA swa; // = wanalysis(tasks[x-1], hps, m, K);
      string fout("results/"+fname+"-res-"+to_string(m)+"-"+to_string(K));
      std::cout << "------> " << fout << "\n";
      ofstream fouts;
      fouts.open(fout, ios::app);
      fouts << "hps size: " << hps.size() << ", (m, K) " << m << ", " << K << ", bounded = " << swa.bounded << ", time = " << swa.time << ", R=" << tasks[x-1].wcrt << ", bp=" << tasks[x-1].bp << endl;
      fouts << "----------*--------------" << endl;
      for ( int j = 1; j <= x; j++)
        fouts << tasks[j-1].get_wcet() << ", " << tasks[j-1].get_period() << ", " << tasks[j-1].get_wcrt() << ", " << tasks[j-1].get_bcrt() << endl;
      fouts << "----------*-------------------";
      fouts.close();

      // to print out bounded/un-bounded results separately
      fouts.open(fout+"-"+to_string(swa.bounded), ios::app);
      fouts << swa.time << std::endl;
      fouts.close();

      if( nsched == 1000) return;
    }

  }

}

//void simul(const int n, const int m, const int nSets, const string& fname) {
//
//  ifstream input(fname);
//  int count = 0;
//  int csched1 = 0, csched2 = 0;
//
//  while (count ++ < nSets) {
//    cout << count << endl;
//    vector<Task> tasks;
//    string line;
//
//    for ( int i = 0; i < n; i ++) {
//      getline(input, line); // read a line from the input file
//      stringstream linestream(line);
//      string data;
//      getline(linestream, data, '\n');
//
//      vector<string> tokens = split(data, ' ');
//
//      double u;
//      int c, p;
//      u = atof(tokens[0].c_str());
//      c = atof(tokens[1].c_str());
//      p = atof(tokens[2].c_str());
//
//      default_random_engine generator(seed++);
//      uniform_real_distribution<double> disD(0.9*p, p);
//      int d = disD(generator);
//      if ( d < c ) d = p;
//
//      Task t(c, d, p);
//      
//      tasks.push_back(t);
//    }
//    getline(input, line); // between two tasksets, there is a line of delimiter
//
//    sort(tasks.begin(), tasks.end(), 
//        [](const Task& a, const Task& b) {return a.get_dline() < b.get_dline();}
//        );
//
//    bool sched_rta = true, bl = true;
//    sched_rta = RTA_LC(tasks, m);
//    if (not sched_rta)
//      bl = baseline(tasks,m);
//
//    if (sched_rta) csched1 ++;
//    if (bl) csched2 ++;
//
//    // schedulability analysis
//    string fout(fname); fout.append(".res.dat");
//    ofstream fouts;
//    fouts.open(fout, ios::app);
//    fouts << count << "  " << csched1 << "  " << csched2 << endl;
//    //fouts << "-------------------------" << endl;
//    //for ( int j = 1; j <= x; j++)
//    //  fouts << tasks[j-1].get_wcet() << ", " << tasks[j-1].get_period() << endl;
//    //fouts << "------------------------------";
//    fouts.close();
//
//
//
//  }
//
//}
//
//
//
///**************** Incrementally generating the taskset ************/
//Task gen_a_task(double Umin=0.05, double Umax=0.45) {
//
//  int Tmin = 1, Tmax = 200;
//
//  int p = 0, c = 0, d = 0;
//  while ( true) {
//    srand (seed ++);
//    p = rand() % (Tmax - Tmin) + Tmin;
//
//    double dc = Umin;
//    if (Umin != Umax) {
//      default_random_engine generator(seed++);
//      uniform_real_distribution<double> disC(Umin, Umax);
//      dc = disC(generator);
//      //if (dc <= Umin) dc = Umin; //continue; 
//      //if (dc > Umax ) dc = Umax; //continue; 
//    }
//    c = dc * p;
//
//    d = p;
//    if ( c > 0 and c < p and c <= d) break;
//  }
//
//  Task t(c, d, p);
//
//  return t;
//
//}

