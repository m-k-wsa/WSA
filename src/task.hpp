#ifndef __TASK_HPP__
#define __TASK_HPP__

#include <string>
#include <vector>

#include "has_unique_id.hpp"


class Task : public HasUniqueId<Task> {

  private:
    std::string name;

  public:

    double wcet, bcet, dline, period;
    double wcrt, bcrt;
    double bp;
    double end_r;
    double jitter;
    double period_lb, period_ub;


    Task ();
    Task (const double c, const double d, const double p); 
    Task (const double c_lb, const double c_ub, const double d, const double p);
    Task (const Task& t);
    Task (const double jitter, const double c_lb, const double c_ub, const double d, const double p_min, const double p_max);


    void set_wcet(const double w) { wcet = w;}
    void set_wcrt(const double w) { wcrt = w;}
    void set_bcrt(const double b) { bcrt = b;}
    void set_bp(const double b) { bp = b;}
    void set_end_r(const double b) { end_r = b;}

    double get_wcet() const { return wcet;}
    double get_bcet() const { return bcet;}
    double get_dline() const { return dline;}
    double get_period() const { return period;}
    double get_wcrt() const { return wcrt;}
    double get_bcrt() const { return bcrt;}
    double get_bp() const { return bp;}
    double get_end_r() const { return end_r;}

    void print() const ;

};


void RTA( Task& tk, const std::vector<Task>& hps);
void compute_bp(Task& tk, const std::vector<Task>& hps);

int NCW(const Task& ti, const int L);
int minWorkload(const Task& ti, const int L);

void compute_wcrt(Task& tk, const std::vector<Task>& hps);
void compute_bcrt(Task& tk, const std::vector<Task>& hps);

double computeIdle(const std::vector<Task> & hps, const int L);
bool BRTA(Task& tk, const std::vector<Task>& hps);

double BI(const double y, const std::vector<Task> &hps);
void compute_bcrt2(Task& tk, const std::vector<Task>& hps);

void tasks_sorting(std::vector<Task> &tks);

std::vector<Task> read_a_taskset(const std::string &fname, const int n);

#endif
