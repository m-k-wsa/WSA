#include "weak.hpp"
#include <iostream>
#include <cmath>
#include <utility>
#include <glpk.h>



using namespace std;

namespace Scan {

  void hello_glpk() {
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
  }







int RTA(const Task& tk, const vector<Task>& hps) {

  int c = tk.get_wcet(), p = tk.get_period();

  int L = c;

  while ( true) {
    int I = ceil(L*1.0 / p) * c;
    for ( auto & x : hps) 
      I += ceil(L*1.0 / x.get_period()) * x.get_wcet();

    if ( I == L) return L;
    else L = I;

  }

}

int compute_idle_t(const std::vector<Task>& hps, const std::vector<pair<int, int> >& idles) {

  int idle_t = 0, w = 0;
  for ( auto & x : idles)
    idle_t += (x.second - x.first);

  for ( auto & x : hps) {

    int wi = idle_t;
    int p = x.get_period(), c = x.get_wcet();

    for ( int o = 0; o <= p; o++) {

      int wio = 0;
      for ( auto & y : idles) {
        int n = (y.second - o) / p;
        int N = ceil( (y.first-o) * 1.0 /p );
        int wii = (n-N) * c;
        if ( wii < 0) wii = 0;
        wio += wii;
      }
      if ( wio < wi) wi = wio;
      if ( wio == 0) break;

    }

    w += wi;
  }

  idle_t -= w;
  if ( idle_t < 0) idle_t = 0;
  return idle_t;

}

int compute_idle_t(const std::vector<Task>& hps, const int start, const int end) {
  int W = 0;
  for ( auto & x : hps) {
    
    int p = x.get_period(), c = x.get_wcet(), d = x.get_dline();
    int w = end - start;
    
    for ( int o = 0; o <= p; o ++)  {
      if (o > start) break;
      int w1 = 0;
      int s = o;
      while (s < end) {
        if ( s >= start and s + d <= end) w1 += c;
        s += p;
      }
      if ( w1 < w) w = w1;
      if ( w1 == 0) break;

    }
    W += w;
  }
  return fmax(0, end-start-W);
}

int wa(const Task& tk, const int K, const std::vector<Task>& hps) {

  int dmiss = 0;
  vector<pair<int, int> > idles;
  int p = tk.get_period(), d = tk.get_dline(), c = tk.get_wcet();

  int idle_t = 0;
  for ( int k = 1; k <= K; k++) {

    idle_t = compute_idle_t(hps, idles);
    int start = (k-1) * p, dline = start + d;

    int L = k * c;

    while (true) {

      int I = k * c + idle_t;
      for ( auto & x : hps) 
        I += ceil(L*1.0 / x.get_period()) * x.get_wcet();
      
      if ( I == L) break;
      else L = I;

    }

    cout << "idle_t = " << idle_t << endl;
    cout << "dline - L = " << dline - L << endl;
    if ( L > dline) dmiss ++;
    //if ( L < start + p) idles.push_back(make_pair(L, start+p)); 
    idles.push_back(make_pair(start+c, start+p)); 
    //if ( L < start + p) idle_t += compute_idle_t(hps, L, start+p);
  }

  return dmiss;
}


}
