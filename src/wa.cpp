#include "wa.hpp"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <iostream>
#include <string>
#include <limits>
using namespace std;

std::string  convert_to_string(const int Number) {
  stringstream convert; 

  convert << Number;

  string Result = convert.str();
  return Result;
}

/***
 * The weakly hard real-time analysis for periodic task systems.
 * Give 'K' successive activations of "ti", can the number of deadline misses exceed 'm'?
 */
statWA wanalysis(const Task& ti, const std::vector<Task>& hps, const int m, const int K) {

  cout << "Weakly hard real-time analysis of a task ti" << endl;
  double ci = ti.get_wcet(), di = ti.get_dline(), pi_min = ti.period_lb, pi_max=ti.period_ub, Ri = ti.get_wcrt(), bri=ti.get_bcrt();
  double BP = ti.get_bp(); int Ni = ceil(BP / pi_max);
  //double ui = ci / pi;
  
  // the vector of wcet, dline, period and utilization of each higher priority task
  vector<double> c, d, p, u, br, jitter;
  for ( int j = 1; j <= hps.size(); j++) {
    double cj = hps[j-1].get_wcet(), dj = hps[j-1].get_dline(), pj = hps[j-1].get_period();
    double jj=hps[j-1].jitter;
    double brj=hps[j-1].get_bcrt();
    c.push_back(cj); d.push_back(dj); p.push_back(pj); u.push_back(cj/pj);
    br.push_back(brj);
    jitter.push_back(jj);
  }
  for (auto &jj : jitter)
    std::cout << "jitter: " << jj << "\n";
  //double totU = ui;
  //for ( int j = 1; j <= hps.size(); j++)
  //  totU += c[j-1]/p[j-1];
  const double M = 10000000, sigma = 0.001;

  vector<Task> copy_hps(hps);
  copy_hps.push_back(ti);
  vector<double> idle_lb; // minimum level-i idle time w.r.c time lengths of pi, 2*pi, 3*pi, ..., K*pi
  for ( int i = 1; i <= K; i++) {
    idle_lb.push_back( fmax(0, computeIdle(copy_hps, i*pi_min)));
    cout << idle_lb[i-1] << endl;
  }



  IloEnv env;

  try {

    IloModel model(env);

    /**
     * ASSUMPTIONS
     *
     * A1: the job in prior to the 1st job is schedulable (otherwise, we can always advance the job windows without reducing
     *      the number of deadline misses)
     *
     * A2: the 1st job must misses its deadline, otherwise, we can always postpone the beginning of the problem windows without
     *      reducing the number of deadline misses)
     *
     */


    /*********************** VARIABLES **********************************/

    // the period is a constant in between 'pi_min' and 'pi_max'
    IloNumVar pi(env, pi_min, pi_max, ILOFLOAT, "period");

    // (a) Busy windows <C>
    IloNumVarArray L(env, 0, 0, IloInfinity, ILOFLOAT);
    for ( int k = 1; k <= K ;k++) {
      string name_bp = "bp" + convert_to_string(k);
      L.add(IloNumVar(env, 0, pi_max-bri, ILOFLOAT, name_bp.c_str())); // [1]
      model.add(IloRange(env, -IloInfinity, L[k-1]-pi, -bri)); // **jitter
    }
    // (b) offsets <C> ***
    IloNumVarArray alpha(env, 0, 0, IloInfinity, ILOFLOAT);
    for ( int j = 1; j <= hps.size() ;j++) {
      string name_alpha = "alpha" + convert_to_string(j);
      alpha.add(IloNumVar(env, 0, p[j-1]-br[j-1]+jitter[j-1], ILOFLOAT, name_alpha.c_str())); //[2]
    }

    // (c) finish times <C>
    IloNumVarArray f(env, 0, 0, IloInfinity, ILOFLOAT);
    for ( int k = 1; k <= K ;k++) {
      string name_ft = "ft" + convert_to_string(k);
      f.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_ft.c_str()));
      
      IloExpr rk(env); rk += L[0] + (k-1)*pi;
      model.add(IloRange(env, bri, f[k-1]-rk, Ri)); // [3]
      if (k > 1)
      {
        //model.add(IloRange(env, ci, f[k-1]-f[k-2], Ri+pi-bri)); // [4]
        model.add(IloRange(env, ci, f[k-1]-f[k-2], IloInfinity)); // [4]
        model.add(IloRange(env, -IloInfinity, f[k-1]-f[k-2]-pi, Ri-bri)); // [4]
      }
    }

    // (d) idle times <C>
    IloNumVarArray idle(env, 0, 0, IloInfinity, ILOFLOAT);
    for ( int k = 1; k <= K ;k++) {
      string name_idle = "idle" + convert_to_string(k);
      idle.add(IloNumVar(env, 0, pi_max-bri, ILOFLOAT, name_idle.c_str())); // [5]
      model.add(IloRange(env, -IloInfinity, idle[k-1]-pi, -bri)); // ***jitter
    }

    // (e) schedulability <B> 
    // b == 1 --> the job is not schedulable; b == 0 --> schedulable
    IloNumVarArray b(env, 0, 0, 1, ILOBOOL);
    for ( int k = 1; k <= K; k++) {
      string name = "b" + convert_to_string(k);
      b.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));

      IloExpr dk(env); dk += L[0] + (k-1)*pi + pi; //di;
      model.add (IloRange(env, 0, M*b[k-1] + dk - f[k-1], IloInfinity)); // [6.l]
      model.add (IloRange(env, -IloInfinity, dk-f[k-1] + sigma - M*(1-b[k-1]), 0)); // [6.r]
    }
    // b[0] == 1: J_1 misses its deadline
    model.add(IloRange(env, 1, b[0], 1)); // [7]

    // (f) interference from the prior job <B>
    //IloNumVarArray bb(env, 0, 0, 1, ILOBOOL);
    for ( int k = 1; k < K; k++) {
      //string name = "bb" + convert_to_string(k);
      //bb.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));

      //IloExpr rkp1(env); rkp1 += L[0] + k*pi;
      //model.add(IloRange( env, 0, rkp1 - f[k-1] + M*bb[k-1], IloInfinity)); // [8.l]
      //model.add(IloRange( env, -IloInfinity, rkp1 - f[k-1] + sigma - M*(1-bb[k-1]), 0)); // [8.r]
      ////idle > 0 --> fbb == 0
      //model.add ( IloRange(env, 0, -idle[k-1] + M*(1-bb[k-1], IloInfinity))); // [9]
      model.add ( IloRange(env, 0, -idle[k-1] + M*(1-b[k-1], IloInfinity))); // [9]
    }

    // (g) counting jobs from a higher priority task <C>
    //// nf_j_k: #jobs of tau_j in [0,f_k)
    //// nL_j_k: #jobs of tau_j in: [0, r_k-L_k) if b_{k-1}=0; [0,f_{k-1}) if  b_{k-1}=1
    std::vector<IloNumVarArray> nf, nL;
    for ( int j = 1; j <= hps.size(); j++) {
      nf.push_back( IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
      nL.push_back( IloNumVarArray(env, 0, 0, IloInfinity, ILOFLOAT));
      for ( int k = 1; k <= K; k++) {
        string name_f = "nf_" + convert_to_string(j) + "_" + convert_to_string(k);
        nf[j-1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_f.c_str()));

        string name_b = "nL_" + convert_to_string(j) + "_" + convert_to_string(k);
        nL[j-1].add(IloNumVar(env, 0, IloInfinity, ILOFLOAT, name_b.c_str()));

      }
    }

    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 1; k <= K; k++) {
        //// a coarse upper bound (uf) [Ex]
        model.add ( IloRange(env,     -IloInfinity,     nf[j-1][k-1] - ceil( (pi_max-bri+(k-1)*pi_max+Ri+jitter[j-1])/p[j-1]),             0) ); //*** jitter

        //// nf [10]
        model.add ( IloRange(env,     0,                nf[j-1][k-1] - (f[k-1] - alpha[j-1]-jitter[j-1])/p[j-1],              IloInfinity) );
        model.add ( IloRange(env,     -IloInfinity,     nf[j-1][k-1] + sigma - (f[k-1] - alpha[j-1]+jitter[j-1])/p[j-1] - 1,  0) ); //****jitter


        if (k > 1) {
          //// a coarse upper bound (uf) [Ex]
          //model.add ( IloRange(env,     -IloInfinity,     nL[j-1][k-1] - ceil( (pi-bri+(k-1)*pi)/p[j-1]) - bb[k-2]*M,             0) );
          model.add ( IloRange(env,     -IloInfinity,     nL[j-1][k-1] - ceil( (pi_max-bri+(k-1)*pi_max+jitter[j-1])/p[j-1]) - b[k-2]*M,             0) ); //***jitter

          //// nL [11]
          //model.add ( IloRange(env,     0,             bb[k-2]*M + nL[j-1][k-1] - (L[0]+(k-1)*pi - L[k-1] - alpha[j-1])/p[j-1],   IloInfinity) );
          //model.add ( IloRange(env,     -IloInfinity,  -bb[k-2]*M + nL[j-1][k-1] + sigma - (L[0]+(k-1)*pi - L[k-1] - alpha[j-1])/p[j-1] - 1,             0) );
          model.add ( IloRange(env,     0,             b[k-2]*M + nL[j-1][k-1] - (L[0]+(k-1)*pi - L[k-1] - alpha[j-1]-jitter[j-1])/p[j-1],   IloInfinity) );
          model.add ( IloRange(env,     -IloInfinity,  -b[k-2]*M + nL[j-1][k-1] + sigma - (L[0]+(k-1)*pi - L[k-1] - alpha[j-1]+jitter[j-1])/p[j-1] - 1,             0) ); //***jitter
        }
        else { // [12]
          model.add ( IloRange(env,     0,     nL[j-1][k-1],   0) );
        }

      }
    }

    // [13]
    for ( int k = 2; k <= K; k++) {
      for ( int j = 1; j <= hps.size(); j++) {
        //model.add(IloRange(env, 0, nL[j-1][k-1] - nf[j-1][k-2] + (1-bb[k-2])*M, IloInfinity));
        //model.add(IloRange(env, -IloInfinity, nL[j-1][k-1] - nf[j-1][k-2] - (1-bb[k-2])*M, 0));
        // b_{k-1}==1 --> nf_{k-1}=nL_k
        model.add(IloRange(env, 0, nL[j-1][k-1] - nf[j-1][k-2] + (1-b[k-2])*M, IloInfinity));
        model.add(IloRange(env, -IloInfinity, nL[j-1][k-1] - nf[j-1][k-2] - (1-b[k-2])*M, 0));
      }
    }



    // (h) refining the job counting for a higher priority task
    
    //// nfb_j_k[]
    vector< std::vector<IloNumVarArray> > nfb;
    //double ub=fmin( Ri, BP-ci);
    //double ubf=Ri;
    double ubf=Ri+pi_max-bri;//***jitter
    //double ub=fmin( Ri+pi-ci, BP-ci);
    for ( int j = 1; j <= hps.size(); j++) {
      vector<IloNumVarArray> v;
      for ( int k = 1; k <= K; k++) {
        IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
        //for ( int i = 1; i <= ceil(Ri/p[j-1]); i++) {
        for ( int i = 1; i <= ceil((ubf+jitter[j-1])/p[j-1]); i++) {
          string name = "nfb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
          vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
        }
        v.push_back( vv);
      }
      nfb.push_back(v);
    }
    ////// nL_j_k + Delta = nf_j_k [14]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 1; k <= K; k++) {
        IloExpr Delta(env);
        //for ( int i = 1; i <= ceil(Ri/p[j-1]); i++)
        for ( int i = 1; i <= ceil((ubf+jitter[j-1])/p[j-1]); i++)
          Delta += nfb[j-1][k-1][i-1];
        model.add(IloRange( env, 0, nL[j-1][k-1] + Delta - nf[j-1][k-1], 0));
      }
    }
    ////// precedent constraint on ""nfb"s [15]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 1; k <= K; k++) {
        //for ( int i = 1; i <= ceil(Ri/p[j-1]); i++) {
        for ( int i = 1; i <= ceil((ubf+jitter[j-1])/p[j-1]); i++) {
          if ( i > 1)
            model.add( IloRange(env, 0, nfb[j-1][k-1][i-1]-nfb[j-1][k-1][i-2], 1));
        }
      }
    }

    //// nLb_j_k[]
    vector< std::vector<IloNumVarArray> > nLb;
    for ( int j = 1; j <= hps.size(); j++) {
      vector<IloNumVarArray> v;
      for ( int k = 1; k <= K; k++) {
        IloNumVarArray vv(env, 0, 0, 1, ILOBOOL);
        for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++) { //***jitter
          string name = "nLb_" + convert_to_string(k) + "_" + convert_to_string(j) + "_" + convert_to_string(i);
          vv.add(IloNumVar(env, 0, 1, ILOBOOL, name.c_str()));
        }
        v.push_back( vv);
      }
      nLb.push_back(v);
    }
    ////// if bb[k-2] == 0, nf[k-2]+Delta_p == nL[k-1] // [16]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 2; k <= K; k++) {
        IloExpr Delta_p(env);
        for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++) //**jitter
          Delta_p += nLb[j-1][k-1][i-1];
        model.add(IloRange( env, 0, nf[j-1][k-2] + Delta_p - nL[j-1][k-1] + b[k-2]*M, IloInfinity));
        model.add(IloRange( env, -IloInfinity, nf[j-1][k-2] + Delta_p - nL[j-1][k-1] - b[k-2]*M, 0));
      }
    }
    ////// if bb[k-2] == 1, nLb == 0 // [17]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 2; k <= K; k++) {
        for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++) { //***jitter
          model.add(IloRange( env, 0, nLb[j-1][k-1][i-1] + (1-b[k-2])*M, IloInfinity));
          model.add(IloRange( env, -IloInfinity, nLb[j-1][k-1][i-1] - (1-b[k-2])*M, 0));
        }
      }
    }
    ////// nLb[j][k][i] >= nLb[j][k][i-1] // [18]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int k = 1; k <= K; k++) {
        for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++) { //***jitter
          if ( i > 1)
            model.add( IloRange(env, 0, nLb[j-1][k-1][i-1]-nLb[j-1][k-1][i-2], 1));
        }
      }
    }
    ////// the cumulative sum of boolean indications [19]
    for ( int j = 1; j <= hps.size(); j++) {
      IloExpr totN(env);
      for ( int k = 1; k <= K; k++) {

        if ( k > 1) {
          for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++)//***jitter
            totN += nLb[j-1][k-1][i-1];
        }

        //for ( int i = 1; i <= ceil(Ri/p[j-1]); i++)
        for ( int i = 1; i <= ceil((ubf+jitter[j-1])/p[j-1]); i++)
          totN += nfb[j-1][k-1][i-1];
        
        model.add( IloRange(env, 0, totN - nf[j-1][k-1], 0));

      }
    }
    ////// nLb_j_1_p = 0 [Ex]
    for ( int j = 1; j <= hps.size(); j++) {
      for ( int i = 1; i <= ceil((pi_max-bri+jitter[j-1])/p[j-1]); i++) { //**jitter
        model.add( IloRange(env, 0, nLb[j-1][0][i-1], 0));
      }
    }

    /*********************** The end of declaration of VARIABLES **********************************/



    /*********************** Weakly hard real-time schedulability analysis **********************************/
    // total number of deadline misses
    IloExpr nDmiss(env);
    for ( int k = 1; k <= K; k++)
      nDmiss += b[k-1];
    model.add(IloRange(env, m+1, nDmiss, K));
    //model.add(IloMaximize(env, nDmiss));
    /********************************************************************************************************/



    /*********************** CONSTRAINTS on schedulability analysis **********************************/




    // C1: constraints on the minimum idle time/workload
    for ( int ni = 1; ni <= K; ni ++) {
      for ( int k = 1; k <= K - ni + 1; k++ ){
        IloExpr idles(env);

        for ( int i = k; i <= k+ni-1; i++) {
            idles += idle[i-1];
        }

        model.add( IloRange(env, idle_lb[ni-1], idles, IloInfinity));

      }
    }


    // C2: idle time inside a job window
    for ( int k = 2; k<= K; k++) {
      IloExpr work(env);
      for ( int j = 1; j <= hps.size(); j++) {
        work += (nL[j-1][k-1] - nf[j-1][k-2]) * c[j-1];
      }
      IloExpr resp(env);
      resp += f[k-2]-L[0]-(k-2)*pi;
      model.add (IloRange(env, -IloInfinity, work + idle[k-2] - b[k-2]*M - (pi - L[k-1] - resp),           0));
      model.add (IloRange(env,            0, work + idle[k-2] + b[k-2]*M - (pi - L[k-1] - resp), IloInfinity));
    }


    /*** the formulation of "finish time", by its precedent job ***/
    /**
     * C4: bb_{k-1} == 1
     **/
    for ( int k = 2; k <= K; k++) {
      IloExpr totI(env);
      for ( int j = 1; j <= hps.size(); j++) {
        totI += (nf[j-1][k-1] - nf[j-1][k-2]) * c[j-1];
      }
      totI += ci;

      model.add( IloRange(env, -IloInfinity, totI - (f[k-1]-f[k-2]) - M*(1-b[k-2]), 0));
      model.add( IloRange(env, 0, totI - (f[k-1]-f[k-2]) + M*(1-b[k-2]), IloInfinity));

    }
    /**
     * C3: b_{k-1} == 0
     **/
    for ( int k = 1; k <= K; k++) {
      IloExpr totI(env);
      for ( int j = 1; j <= hps.size(); j++) {
        totI += (nf[j-1][k-1] - nL[j-1][k-1]) * c[j-1];
      }
      totI += ci;

      IloExpr refer(env); refer = L[0]+(k-1)*pi-L[k-1];
      if ( k == 1) {
        model.add( IloRange(env, -IloInfinity, totI - (f[k-1]-refer), 0));
        model.add( IloRange(env, 0, totI - (f[k-1]-refer), IloInfinity));
      }
      else {
        model.add( IloRange(env, -IloInfinity, totI - (f[k-1]-refer) - M*(b[k-2]), 0));
        model.add( IloRange(env, 0, totI - (f[k-1]-refer) + M*(b[k-2]), IloInfinity));
      }

    }

    /**
     * C5: "ft" by accummulating previous workload and idles
     **/
    for ( int k = 1; k <= K; k++) {
      IloExpr totI(env);
      for ( int j = 1; j <= hps.size(); j++) {
        totI += nf[j-1][k-1] * c[j-1];
      }
      totI += k*ci;

      IloExpr totIdle(env);
      for ( int i = 1; i < k; i++) totIdle += idle[i-1];

      model.add( IloRange(env, 0, totI + totIdle - f[k-1], 0));
    }

    /************************** C.12 *************************/
    // C6: let's try to refine the beginning of "bp"
    for ( int k = 2; k <= K; k++) {
      for ( int j = 1; j <= hps.size(); j ++) 
        model.add(IloRange(env, -IloInfinity, alpha[j-1] + (nL[j-1][k-1]-1)*p[j-1] + br[j-1] + sigma - (L[0]+(k-1)*pi-L[k-1])-b[k-2]*M, 0));
    }
    // C7: let's try to refine the beginning of "ft"
    for ( int k = 1; k <= K; k++) {
      for ( int j = 1; j <= hps.size(); j++) {
        model.add(IloRange(env, -IloInfinity, alpha[j-1] + (nf[j-1][k-1]-1)*p[j-1] + br[j-1] + sigma - f[k-1], 0));
      }
    }


    // C8: the longest bust period

    for ( int k = 1; k <= K - Ni + 1; k++ ){
      IloExpr sumE(env);

      //if ( k < K)
        sumE +=  L[k-1]+pi  - (1-b[k-1])*M;
      //else 
      //  sumE +=  L[k-1]+pi;

      for ( int i = k+1; i <= k+Ni-1; i++) {
        if ( i < k+Ni-1)
          sumE += pi - (1-b[i-1])*M;
        else
          sumE += f[i-1] - (L[0]+(i-1)*pi);
      }
      model.add( IloRange(env, -IloInfinity, sumE, BP));
    }






    //// Extracting the model
    IloCplex cplex(env);

    const int timeLimit = 60*30;
    //const int timeLimit = 7200; //60*30*4;
    cplex.setParam(IloCplex::TiLim,timeLimit);
    //cplex.setParam(IloCplex::TreLim,1024);


    double ss = cplex.getCplexTime();



    cplex.extract(model);
    bool fea = cplex.solve();

    double ee = cplex.getCplexTime();
    bool bounded = false;
    statWA swa;
    swa.time = ee - ss;
  

    if (fea) {
      swa.bounded = false; 
    }
    else {
      if (swa.time < timeLimit)
        swa.bounded = true;
      else
        swa.bounded = false;
    }
    
    env.end();
    return swa;
  
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }

  env.end();


}
