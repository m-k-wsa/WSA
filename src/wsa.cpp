#include <unistd.h>
#include "task.hpp"
#include "wa.hpp"
#include "exp.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>

#include <iostream>
using namespace std;

/*funcion that show the help information*/
void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-n  number of tasks"<<endl;
  cout<<"         "<<"-m  the m parameter"<<endl;
  cout<<"         "<<"-k  the K parameter"<<endl;
  cout<<"         "<<"-f  the file containing input tasks in the form of (C, D, T)"<<endl;
  cout<<"         "<<"-x  kill the deadline-miss job"<<endl;
}


int main (int argc,char *argv[])
{
  int n=-1, m=-1, K=-1;
  string fname="";
  bool job_kill=false;

  char tmp;
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if(argc == 1)
  {
    showhelpinfo(argv[0]);
    exit(1);
  }
  while((tmp=getopt(argc,argv,"xh:n:m:k:f:"))!=-1)
  {
    switch(tmp)
    {
      /*option h show the help infomation*/
      case 'h':
        showhelpinfo(argv[0]);
        break;
      case 'n':
        n=atoi(optarg);
        continue;
      case 'm':
        m=atoi(optarg);
        continue;
      case 'k':
        K=atoi(optarg);
        continue;
      case 'x':
        job_kill=true;
        continue;
      case 'f':
        fname=string(optarg);
        break;
      default:
        showhelpinfo(argv[0]);
        break;
    }
  }
  if(n<0 or m<0 or K<0 or fname=="")
  {
    showhelpinfo(argv[0]);
    exit(1);
  }

  vector<Task> tasks=read_a_taskset(fname, n);

  for ( int x = 1; x <= n; x++) {
    vector<Task> hps;
    for ( int y = 1; y < x; y++) 
      hps.push_back( tasks[y-1]);
    compute_wcrt(tasks[x-1], hps);
    compute_bp(tasks[x-1], hps);
    compute_bcrt(tasks[x-1], hps);

    tasks[x-1].print();
    if(tasks[x-1].bcrt>tasks[x-1].dline) continue;
    if(tasks[x-1].wcrt<=tasks[x-1].dline) continue;
    std::cout << "+++++++++++++++++++++++++++++++\n";
    std::cout << "Task " << tasks[x-1].get_id() << " is un-schedulable\n";
    std::cout << "Running weakly-hard analysis: m=" << m << ", K=" << K << " ... " << std::endl;
    std::cout << "job_kill: " << job_kill << endl;

    // streambuf* orig_buf = cout.rdbuf();
    // cout.rdbuf(NULL);

    statWA swa;
    if(job_kill)
      swa=wanalysis_kill(tasks[x-1], hps, m, K);
    else
      swa = wanalysis(tasks[x-1], hps, m, K);

    // cout.rdbuf(orig_buf);

    if(swa.bounded>0) std::cout << "The (m, K) property holds\n";
    else std::cout << "The (m, K) property does not hold\n";
    std::cout << "Time = " << swa.time << " sec " << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++\n";
  }

  return 0;
}
