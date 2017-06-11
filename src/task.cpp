#include "task.hpp"
#include "util.hpp"

#include <iostream>
#include <limits>
#include <cmath>
#include <utility>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

using namespace std;

void tasks_sorting(std::vector<Task> &tks)
{
  std::cout << "tasks sorting\n";
  int n = tks.size();
  for ( int c=0; c < n-1; c++)
  {
    cout << "c: " << c << std::endl;
    for(int d=0; d < n-c-1; d++)
    {
      cout << "-->d: " << d << std::endl;
      if(tks[d].get_dline() > tks[d+1].get_dline())
      {
        cout << "-->--> before swap" << std::endl;
        Task tmp = tks[d];
        tks[d] = tks[d+1];
        tks[d+1] = tmp;
        cout << "-->--> after swap" << std::endl;
      }
      else cout << "-->--> no swap\n";
    }
  }
  std::cout << "what happens here\n";
}

template<>
int HasUniqueId<Task>::counter = 0;

Task::Task() : wcet(0), dline (0), period (0), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}

Task::Task (const double c, const double d, const double p)
  : wcet(c), dline(d), period(p), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}

Task::Task (const double c_lb, const double c_ub, const double d, const double p)
  : bcet(c_lb), wcet(c_ub), dline(d), period(p), wcrt(numeric_limits<double>::max()), bcrt(numeric_limits<double>::min()) {}


Task::Task (const Task& t) 
  : HasUniqueId<Task>(t), wcet(t.wcet), dline(t.dline), period(t.period), wcrt(t.wcrt), bcrt(t.bcrt), bcet(t.bcet) {}

void Task::print() const {
  cout << "task id: " << get_id() << ", (" << bcet << ", " << wcet << ", " << get_dline() << ", " << period << ") ";
  cout << "wcrt = " << wcrt << ", bcrt = " << bcrt;
  cout << endl;
}

void RTA(Task& tk, const vector<Task>& hps) {

  int c = tk.get_wcet(), p = tk.get_period();

  int bp = 0;
  int L = c;
  int wcrt = 0;
  int k = 1; bool flag = true;
  while ( flag) {
    int X = k * c;

    while ( true) {
      int I = k * c;
      for ( int i = 0; i < hps.size(); i++)
        I += ceil(X*1.0/hps[i].get_period()) * hps[i].get_wcet();

      if ( I > X) X = I;
      else {
        if ( I - (k-1)*p > wcrt) wcrt = I - (k-1)*p;
        if (I <= (k-1)*p + p) { bp = I; flag = false;}
        else k++;
        break;
      }

    }

  }
  tk.set_wcrt(wcrt);
  //tk.set_bp(bp);
}

void compute_wcrt(Task& tk, const vector<Task>& hps) {

  double c = tk.get_wcet(), p = tk.get_period();

  double wcrt = 0;
  int k = 0; bool flag = true;
  while ( flag) {
    k++;
    double X = k * c;
    while ( true) {
      double I = k*c;
      for ( int i = 0; i < hps.size(); i++)
      {
        //cout << "I: " << I << std::endl;
        //cout << "X: " << X << endl;
        //cout << "period: " << hps[i].get_period() << endl;
        //cout << "before ceil: " << X*1.0/hps[i].get_period() << endl; 
        //cout << "after ceil: " << ceil(X*1.0/hps[i].get_period()) << endl; 
        //cout << "wcet: " << hps[i].get_wcet() << endl;

        //cout << ceil(X*1.0/hps[i].get_period()) * hps[i].get_wcet() << endl;;
        I += ceil(X*1.0/hps[i].get_period()) * hps[i].get_wcet();
      }

      if ( I > X) X = I;
      else break;

    }
    double tmp=X-(k-1)*p;
    if(tmp>wcrt) wcrt=tmp;
    if(tmp<=p) break;
  }
  tk.set_wcrt(wcrt);
}

void compute_bcrt(Task& tk, const vector<Task>& hps) {

  double c = tk.bcet, p = tk.get_period();
  double X=tk.get_wcrt();
  while(true)
  {
    double I = c;
    for(int i=0; i<hps.size(); i++)
    {
      double p = hps[i].get_period();
      I += ceil( (X-p)/p ) * hps[i].bcet;
    }
    if(I<X) X = I;
    else break;
  }
  tk.set_bcrt(X);
}

void compute_bp(Task& tk, const vector<Task>& hps) {
  double c = tk.get_wcet(), p = tk.get_period();

  double bp = fmax(c, tk.get_wcrt());

  while (true) {
    //std::cout << "bp on growing: " << bp << std::endl;

    double X = 0;
    
    for ( int i = 0; i < hps.size(); i++)
      X += ceil(bp*1.0/hps[i].get_period()) * hps[i].get_wcet();
    X += ceil(bp*1.0/p) * c;

    if( X == bp) {
      break;
      ////double Y = X + 0.1;
      //double Y = X + 0.001;

      //X = ceil(Y*1.0/p) * c;
      //for ( int i = 0; i < hps.size(); i++)
      //  X += ceil(Y*1.0/hps[i].get_period()) * hps[i].get_wcet();

      //if ( X >= Y) bp = X;
      //else break;
    }
    else bp = X;
  }
  tk.set_bp(bp);
}

int NCW(const Task& ti, const int L) {

  int p = ti.get_period(), c = ti.get_wcet();
  int w = (L/p) * c;

  if (L%p <= c) w += L%p;
  else w+= c;

  return w;
}

int minWorkload(const Task& ti, const int L) {

  double p = ti.get_period(), c = ti.get_wcet();
  int w = 0;
  int l = L - (p - c);
  if ( l < 0) return 0;
  else return c/p * l;

  return w;


}

double computeIdle(const std::vector<Task> & hps, const int L) {

  for ( int c = 1; c <= L; c++) {
    Task tk(c, L, L);
    //if ( not BRTA(tk, hps)) return c-1;
    if (BRTA(tk, hps)) continue;
    double delta=0.01;
    double cc=c-1+delta;
    while(cc<=c)
    {
      Task tk2(cc, L, L);
      if(not BRTA(tk2, hps)) return cc-delta;
      cc+=delta;
    }
    return c-1;
  }

}

bool BRTA(Task& tk, const vector<Task>& hps) {

  double c = tk.get_wcet(), p = tk.get_period(), d = tk.get_dline();

  double bp = 0;
  double X = c;
  double wcrt = 0;
  int k = 1; bool flag = true;

  while ( true) {
    double I = c;
    for ( int i = 0; i < hps.size(); i++)
      I += ceil(X*1.0/hps[i].get_period()) * hps[i].get_wcet();

    if ( I > d) return false;
    if ( I > X) X = I;
    else return true;

  }

}

double BI(const double y, const std::vector<Task> &hps)
{
  // comput bu
  double bu=0;
  for ( int i=1; i<=hps.size(); i++)
    bu += hps[i-1].get_bcet()/hps[i-1].get_period();
  // compute wr
  double wr_high=y;
  for (int i=1; i<=hps.size(); i++)
    wr_high+=hps[i-1].get_bcet();

  double wr=wr_high/(1-bu);

  Task t(y, y, wr, wr);

  compute_bcrt(t, hps);
  std::cout << "y = " << y << ", period: " << wr << ", bi: " << t.get_bcrt() << std::endl;
  return t.get_bcrt();
}

void compute_bcrt2(Task& tk, const std::vector<Task>& hps)
{
  int wl=ceil( tk.get_bp()/tk.get_period()); 
  double bcrt=0;
  for ( int i = 1; i <= wl; i++)
  {
    double bi=BI(i*tk.get_bcet(), hps)-(i-1)*tk.get_period(); 
    if(bcrt<bi) bcrt=bi;
  }
  tk.set_bcrt(bcrt);
}

std::vector<Task> read_a_taskset(const string &fname, const int n)
{
  ifstream input(fname);
  vector<Task> tasks;
  string line;

  for ( int i = 0; i < n; i ++) {
    getline(input, line); // read a line from the input file
    stringstream linestream(line);
    string data;
    getline(linestream, data, '\n');
    vector<string> tokens = split(data, ' ');

    double c, d, p;
    c = atof(tokens[0].c_str());
    d = atof(tokens[1].c_str());
    p = atof(tokens[2].c_str());

    Task t(c, c, d, p);
    tasks.push_back(t);
  }

  return tasks;
}

bool com_dline(const Task& a, const Task& b) {return a.get_dline() < b.get_dline();}
