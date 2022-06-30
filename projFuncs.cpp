#include <cmath>
#include <climits>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <numeric>
#include <omp.h>
#include <projFuncs.h>
using namespace std;

// SEE HEADER FILE FOR MORE DOCUMENTATION


// STRUCTS
struct Parameters {
    double eps;
    int maxIter;
    
    vector<double> corr_time;
    double syst_dt;
    double syst_m;
    vector<double> syst_g;
    double syst_umin;
    double syst_umax;
    double syst_vmax;
    double syst_wmax;
    double syst_utheta;
    double corr_vmin;
    double corr_vmax;
    vector<vector<double>> corr_cen;
    vector<vector<double>> corr_dir;
    vector<double> corr_len;
    vector<double> corr_rad;

    vector<double> traj_r0;
    vector<double> traj_rf;
    vector<double> traj_v0;
    vector<double> traj_vf;
    vector<double> traj_wf;
    vector<double> traj_uf;
    double traj_wei_w;

    double omega;
    double rho;
    double lam;
};

struct Trajectory {
    vector<vector<double>> r;
    vector<vector<double>> v;
    vector<vector<double>> u;
    vector<vector<double>> w;
    double solTime;
};

struct minTime {
    double tF;
    double tau;
    double alpha;
    double beta;
    double ratio;
    double eps;
    vector<double> r0;
    vector<double> rf;
    vector<double> vf;
    vector<double> wf;
    vector<double> uf;
    vector<double> cen;
    vector<double> dir;
    double len;
    double rad;
};




// Helper Functions
inline double norm(vector<double> &c) {
    double inner = 0.0;
    for(long unsigned int i = 0; i < c.size(); i++) {
        inner += pow(c[i],2);
    }
    return sqrt(inner);
}

inline void printVector(vector<double> &u) {
    for(long unsigned int i = 0; i < u.size(); i++) {
        cout << u[i] << " ";
    }
    cout << endl;
}

inline void printMatrix(vector<vector<double>> &A) {
    for(long unsigned int i = 0; i < A.size(); i++) {
        for(long unsigned int j = 0; j < A[i].size(); j++) {
            cout << A[i][j] << " ";
        }   
        cout << endl;
    }
}

inline vector<vector<double>> transpose(vector<vector<double>> &A) {
    vector<vector<double>> T(A[0].size(), vector<double>(A.size(),0));
    for(long unsigned int i = 0; i < A.size(); i++) {
        for(long unsigned int j = 0; j < A[i].size(); j++) {
            T[i][j] = A[j][i];
        }
    }
    return T;
}

inline double normPower(vector<vector<double>> &x, vector<vector<double>> &u) {
    double s = 0;
    for(long unsigned int t = 0; t < x.size(); t++) {
        s += pow(norm(x[t]), 2);
    }
    for(long unsigned int t = 0; t < u.size(); t++) {
        s += pow(norm(u[t]), 2);
    }
    return sqrt(s);
}

inline double normDiff(vector<vector<double>> &v, vector<vector<double>> &vBar) {
    double s = 0;
    for(long unsigned int t = 0; t < v.size(); t++) {
        vector<double> temp(v[t].size());
        for(long unsigned int i = 0; i < v[t].size(); i++) {
            temp[i]= v[t][i] - vBar[t][i];
        }
        s += pow(norm(temp), 2);
    }
    return s;
}

inline double powSum(vector<vector<double>> &x, vector<vector<double>> &x1) {
    double sum = 0;
    for(long unsigned int i = 0; i < x.size(); i++) {
        for(long unsigned int j = 0; j < x[i].size(); j++) {
            sum += pow(x[i][j] - x1[i][j], 2);
        }
    }
    return sum;
}

inline double simplePow(vector<vector<double>> &x) {
    double sum = 0;
    for(long unsigned int i = 0; i < x.size(); i++) {
        for(long unsigned int j = 0; j < x[i].size(); j++) {
            sum += pow(x[i][j], 2);
        }
    }
    return sum;
}








// PROJECTION FUNCTIONS
inline vector<double> projPosOrth(vector<double> &c) {
    vector<double> x(c.size());
    for(long unsigned int i = 0; i < c.size(); i++) {
        x[i] = max(c[i], 0.0);
    }
    return x;
}

inline vector<double> projBox(vector<double> &c, vector<double> &l, vector<double> &u) {
    vector<double> x(c.size());
    for(long unsigned int i = 0; i < c.size(); i++) {
        x[i] = min(max(c[i], l[i]), u[i]);
    }
    return x;
}

inline vector<double> projBall(vector<double> &c, double &p) {
    vector<double> x(c.size());
    const double Cnorm = norm(c);
    for(long unsigned int i = 0; i < c.size(); i++) {
        x[i] = (p/max(Cnorm, p)) * c[i];
    }
    return x;
}

inline vector<double> projCone(vector<double> &c, double &g) {
    const int n = c.size();
    const double b = c[n-1];
    vector<double> a(n);
    a = c;
    a.resize(n-1);
    const double Anorm = norm(a);
    vector<double> x(n);

    if(Anorm <= g*b) {
        return c;
    } else if(g*Anorm <= -1.0*b) {
        x.assign(n, 0);
        return x;
    } else {
        const double coeff = (g*Anorm + b)/(g*g + 1.0);
        const double mult = coeff*g/Anorm;
        for(long unsigned int i = 0; i < a.size(); i++) {
            x[i] = mult*a[i];
        }
        x[n-1] = coeff;
        return x;
    }
}

inline vector<double> projConeBall(vector<double> &x, double &theta, double &rho) {
    const double gamma = tan(theta);
    const double beta = x.back();
    vector<double> a(x);
    a.pop_back();
    const double delta = norm(a);

    if(delta > gamma*beta) {
        if(gamma*delta <= -beta) {
            x.assign(x.size(), 0);
        }
        else {
            const double factor = (gamma*delta+beta)/(pow(gamma,2)+1);
            for(int i = 0; i < x.size()-1; i++) {
                x[i] = factor*(gamma/delta)*a[i];
            }
            x.back() = factor;
        }
    }

    const double factor = rho/max(norm(x), rho);
    for(int i = 0; i < x.size(); i++) {
        x[i] = factor * x[i];
    }

    return x;
}

inline vector<double> projCylinder(vector<double> &c, vector<double> &h, 
                                   vector<double> &e, double &eta, double &rho) {
    vector<double> x(c);
    double x1 = 0.0;
    vector<double> x2(c);
    for(long unsigned int i = 0; i < x.size(); i++) {
        x[i] = c[i] - h[i];
        x1 += e[i]*x[i];
    }
    for(long unsigned int i = 0; i < x.size(); i++) {
        x2[i] = x[i] - x1*e[i];
    }
    x1 = max(min(x1, eta), -1.0*eta);
    
    const double factor = rho/max(norm(x2), rho);
    for(long unsigned int i = 0; i < x.size(); i++) {
        x2[i] = factor*x2[i];
        x[i] = h[i] + x1*e[i] + x2[i];
    }

    return x;
}





// LARGER FORM ALGORITHMS
inline double power_iter(const long unsigned int &tau, struct Parameters &param) {
    const double dt = param.syst_dt;
    const double dt2 = param.syst_dt/(2.0*param.syst_m);
    const double dt3 = pow(param.syst_dt,2)/(3.0*param.syst_m);
    const double dt4 = pow(param.syst_dt,2)/(6.0*param.syst_m);

    vector<vector<double>> r(tau+1, vector<double>(3,0));
    vector<vector<double>> v(tau+1, vector<double>(3,0));
    vector<vector<double>> u(tau+1, vector<double>(3,0));
    vector<vector<double>> w(tau+1, vector<double>(3,0));
    vector<vector<double>> cr(tau+2, vector<double>(3,0));
    vector<vector<double>> cv(tau+2, vector<double>(3,0));
    vector<vector<double>> cu(tau+2, vector<double>(3,0));
    vector<double> lu(tau+1);

    for(long unsigned int t = 0; t < r.size(); t++) {
        for(long unsigned int i = 0; i < r[t].size(); i++) {
            r[t][i] = (rand()%10)/100.0;
            v[t][i] = (rand()%10)/100.0;
            u[t][i] = (rand()%10)/100.0;
            w[t][i] = (rand()%10)/100.0;
        }
    }
    
    // sigma initial
    double sigma = pow(simplePow(r) + simplePow(v) + simplePow(u) + simplePow(w), 0.5);
    double sigma1 = 2.0*sigma;
    while(abs(sigma - sigma1) > param.eps) {
        sigma1 = sigma;
        
        for(long unsigned int t = 0; t < r.size(); t++) {
            for(long unsigned int i = 0; i < r[t].size(); i++) {
                r[t][i] = r[t][i]/sigma;
                v[t][i] = v[t][i]/sigma;
                u[t][i] = u[t][i]/sigma;
                w[t][i] = w[t][i]/sigma;
            }
        }
        
        // 1st set of multiplications
        for(long unsigned int t = 0; t < tau; t++) {
            for(long unsigned int i = 0; i < cr[t].size(); i++) {
                cr[t+1][i] = r[t+1][i] - r[t][i] - dt*v[t][i] - dt3*u[t][i] - dt4*u[t+1][i];
                cv[t+1][i] = v[t+1][i] - v[t][i] - dt2*u[t][i] - dt2*u[t+1][i];
                cu[t+1][i] = u[t+1][i] - u[t][i] - w[t][i];
            }
        }

        // 2nd set of multiplications
        for(long unsigned int t = 0; t < tau+1; t++) {
            lu[t] = u[t][2];
            for(long unsigned int i = 0; i < r[t].size(); i++) {
                r[t][i] = cr[t][i] - cr[t+1][i];
                v[t][i] = -dt*cr[t+1][i] + cv[t][i] - cv[t+1][i];
                u[t][i] = -dt4*cr[t][i] - dt3*cr[t+1][i] - dt2*cv[t][i] - dt2*cv[t+1][i] + cu[t][i] - cu[t+1][i] + lu[t];
                w[t][i] = -1.0*cu[t+1][i];
            }
        }

        // compute estimate of Hnorm
        sigma = pow(simplePow(r) + simplePow(v) + simplePow(u) + simplePow(w), 0.5);
    }

    return 1.01*sigma;
}

int pipg_traj(struct Parameters &param, struct minTime &m, struct Trajectory &traj, vector<double> &tau) {
    long unsigned int tF = (int)tau.back();
    
    const vector<double> g = param.syst_g;
    const double dt = param.syst_dt;
    const double dt2 = param.syst_dt/(2.0*param.syst_m);
    const double dt3 = pow(param.syst_dt,2)/(3.0*param.syst_m);
    const double dt4 = pow(param.syst_dt,2)/(6.0*param.syst_m);

    const double sigma = power_iter(tF, param);
    const double alpha = 2.0/(param.lam+pow((param.lam+4.0*param.omega*sigma), 0.5));
    const double beta  = param.omega*alpha;

    vector<vector<double>> r(tF+1, vector<double>(3,0));
    vector<vector<double>> v(tF+1, vector<double>(3,0));
    vector<vector<double>> u(tF+1, vector<double>(3,0));
    vector<vector<double>> w(tF  , vector<double>(3,0));
    for(long unsigned int t = 0; t < r.size(); t++) {
        for(long unsigned int i = 0; i < r[t].size(); i++) {
            r[t][i] = (rand()%600)/100.0-3.0;
            v[t][i] = (rand()%600)/100.0-3.0;
            u[t][i] = (rand()%600)/100.0-3.0;
            if(t < w.size()) {
                w[t][i] = (rand()%600)/100.0-3.0;
            }
        }
    }
    vector<vector<double>> rp(r);
    vector<vector<double>> vp(v);
    vector<vector<double>> up(u);
    vector<vector<double>> wp(w);

    vector<vector<double>> cr(tF+2, vector<double>(3,0));
    vector<vector<double>> cv(tF+2, vector<double>(3,0));
    vector<vector<double>> cu(tF+2, vector<double>(3,0));
    vector<double>         lu(tF+1);
    vector<vector<double>> crp(cr);
    vector<vector<double>> cvp(cv);
    vector<vector<double>> cup(cu);
    vector<double>         lup(lu);

    vector<vector<double>> r1(r);
    vector<vector<double>> v1(v);
    vector<vector<double>> u1(u);
    vector<vector<double>> w1(w);
    vector<vector<double>> cr1(cr);
    vector<vector<double>> cv1(cv);
    vector<vector<double>> cu1(cu);
    vector<double>         lu1(lu);

    const int initIter = 800;
    const int freq = 100;
    int status = 0;
    double ratio = 0;
    double err = 0.0;
    double errOld = 1.0*INT_MAX;
    for(long unsigned int iter = 0; iter < param.maxIter; iter++) {
        cr1 = crp;
        cv1 = cvp;
        cu1 = cup;
        lu1 = lup;

        // primal variable update: gradient+PI
        for(long unsigned int t = 0; t < tF; t++) {
            for(long unsigned int i = 0; i < r[t].size(); i++) {
                rp[t][i] = r[t][i] - alpha*(cr[t][i] - cr[t+1][i]);
                vp[t][i] = v[t][i] - alpha*(-1.0*dt*cr[t+1][i] + cv[t][i] - cv[t+1][i]);
                up[t][i] = u[t][i] - alpha*(u[t][i] - dt4*cr[t][i] - dt3*cr[t+1][i] - dt2*cv[t][i] - dt2*cv[t+1][i] + cu[t][i] - cu[t+1][i] + (i==2)*lu[t]);
                wp[t][i] = w[t][i] - alpha*(param.traj_wei_w*w[t][i] - cu[t+1][i]);
            }
        }
        
        // Initial Conditions
        rp[0] = param.traj_r0;
        vp[0] = param.traj_v0;
        up[0] = projConeBall(up[0], param.syst_utheta, param.syst_umax);
        wp[0] = projBall(wp[0], param.syst_wmax);
        
        // Final Conditions
        rp.back() = param.traj_rf;
        vp.back() = param.traj_vf;
        up.back() = param.traj_uf;
        
        // Projection of Primal Variable
        for(long unsigned int t = 1; t < tF; t++) {
            vp[t] = projBall(vp[t], param.syst_vmax);
            up[t] = projConeBall(up[t], param.syst_utheta, param.syst_umax);
            wp[t] = projBall(wp[t], param.syst_wmax);
        }
        
        for(long unsigned int n = 0; n < tau.size()-1; n++) {
            for(long unsigned int t = tau[n]+1; t < tau[n+1]+1; t++) {
                rp[t] = projCylinder(rp[t], param.corr_cen[n], param.corr_dir[n], param.corr_len[n], param.corr_rad[n]);
            }
        }                

        for(long unsigned int t = 0; t < r1.size(); t++) {
            for(long unsigned int i = 0; i < r1[t].size(); i++) {
                // Compute PI feedback
                r1[t][i] = 2*rp[t][i] - r[t][i];
                v1[t][i] = 2*vp[t][i] - v[t][i];
                u1[t][i] = 2*up[t][i] - u[t][i];

                // Extrapolation for Primal Variables
                r[t][i] = (1-param.rho)*r[t][i] + param.rho*rp[t][i];
                v[t][i] = (1-param.rho)*v[t][i] + param.rho*vp[t][i];
                u[t][i] = (1-param.rho)*u[t][i] + param.rho*up[t][i];

                if(t < w1.size()) {
                    w1[t][i] = 2*wp[t][i] - w[t][i];
                    w[t][i] = (1-param.rho)*w[t][i] + param.rho*wp[t][i];
                }
            }
        }
        
        // Update Dual Variables
        for(long unsigned int t = 0; t < tF; t++) {
            for(long unsigned int i = 0; i < crp[t].size(); i++) {
                crp[t+1][i] = cr[t+1][i] + beta*(r1[t+1][i] - r1[t][i] - dt*v1[t][i] - dt3*u1[t][i] - dt4*u1[t+1][i] - 0.5*pow(dt,2)*g[i]);
                cvp[t+1][i] = cv[t+1][i] + beta*(v1[t+1][i] - v1[t][i] - dt2*u1[t][i] - dt2*u1[t+1][i] - dt*g[i]);
                cup[t+1][i] = cu[t+1][i] + beta*(u1[t+1][i] - u1[t][i] - w1[t][i]);
            }
        }

        // Extrapolation for Dual Variables
        for(long unsigned int t = 0; t < cr.size(); t++) {
            for(long unsigned int i = 0; i < cr[t].size(); i++) {
                cr[t][i] = (1-param.rho)*cr[t][i] + param.rho*crp[t][i];
                cv[t][i] = (1-param.rho)*cv[t][i] + param.rho*cvp[t][i];
                cu[t][i] = (1-param.rho)*cu[t][i] + param.rho*cup[t][i];
            }
        }
        
        err = 0.0;
        for(long unsigned int t = 0; t < tF+1; t++) {
            // Update 'lu' Dual Variable
            lup[t] = min(0.0, lu[t] + beta*(u1[t][2] - param.syst_umin));

            // Extrapolation for 'lu' Dual Variable
            lu[t] = (1-param.rho)*lu[t] + param.rho*lup[t];
            
            // Error Calculation
            err += pow(lup[t] - lu1[t], 2);
        }

        // Check if gotten far enough to estimate termination
        int iter1 = iter - initIter;
        if(iter1 >= 0) {
            if(iter1 == 0) {
                err += powSum(crp,cr1) + powSum(cvp,cv1) + powSum(cup,cu1);
                err = err/(pow(param.rho*beta, 2)*(simplePow(rp)+simplePow(vp)+simplePow(up)+simplePow(wp)));
                err = min(errOld, err);
                errOld = err;

                if(err < m.eps) {
                    status = 1;
                    break;
                }
            }
            else if(iter1 % freq == 0) {
                err += powSum(crp,cr1) + powSum(cvp,cv1) + powSum(cup,cu1);
                err = err/(pow(param.rho*beta, 2)*(simplePow(rp)+simplePow(vp)+simplePow(up)+simplePow(wp)));
                err = min(errOld, err);

                if(err < m.eps) {
                    status = 1;
                    break;
                }
                else {
                    ratio = err/errOld;
                    if(ratio > m.ratio) {
                        status = 0;
                        break;
                    }
                    errOld = err;
                }
            }
        }
    }
    traj.r = rp;
    traj.v = vp;
    traj.u = up;
    traj.w = wp;

    return status;
}

inline double pt2pt(struct Parameters &param, struct minTime &m) {
    const vector<double> g = param.syst_g;
    const long unsigned int maxI = param.maxIter;
    const long unsigned tF = m.tF;
    m.rf = projCylinder(m.rf, m.cen, m.dir, m.len, m.rad);
    m.vf = projBall(m.vf, param.syst_vmax);
    m.uf = projConeBall(m.uf, param.syst_utheta, param.syst_umax);
    m.wf = projBall(m.wf, param.syst_wmax);

    const double dt = param.syst_dt;
    const double dt2 = dt/(2*param.syst_m);
    const double dt3 = pow(dt,2)/(3*param.syst_m);
    const double dt4 = dt3/2;

    vector<vector<double>> r(tF+1, vector<double>(3,0));
    vector<vector<double>> v(tF+1, vector<double>(3,0));
    vector<vector<double>> u(tF+1, vector<double>(3,0));
    vector<vector<double>> w(tF+1, vector<double>(3,0));
    for(long unsigned int t = 0; t < r.size(); t++) {
        for(long unsigned int i = 0; i < r[t].size(); i++) {
            r[t][i] = (rand()%10)/100.0;
            v[t][i] = (rand()%10)/100.0;
            u[t][i] = (rand()%10)/100.0;
            w[t][i] = (rand()%10)/100.0;
        }
    }

    vector<vector<double>> cr(tF+2, vector<double>(3,0));
    vector<vector<double>> cv(tF+2, vector<double>(3,0));
    vector<vector<double>> cu(tF+2, vector<double>(3,0));
    vector<double> lu(tF+1);

    vector<vector<double>> rp(r);
    vector<vector<double>> vp(v);
    vector<vector<double>> up(u);
    vector<vector<double>> wp(w);
    vector<vector<double>> crp(cr);
    vector<vector<double>> cvp(cv);
    vector<vector<double>> cup(cu);
    vector<double> lup(lu);

    vector<vector<double>> r1(rp);
    vector<vector<double>> v1(vp);
    vector<vector<double>> u1(up);
    vector<vector<double>> w1(wp);
    vector<vector<double>> cr1(crp);
    vector<vector<double>> cv1(cvp) ;
    vector<vector<double>> cu1(cup);
    vector<double> lu1(lup);
    int initIter = 1200;
    int freq = 150;
    double err = INT_MAX;
    double errOld;
    double ratio = -1;
    
    for(long unsigned int iter = 0; iter < maxI; iter++) {
        // Primal Variable Update: Gradient+PI
        for(long unsigned int t = 0; t < r.size(); t++) {
            for(long unsigned int i = 0; i < r[t].size(); i++) {
                rp[t][i] = r[t][i] - m.alpha*(cr[t][i] - cr[t+1][i]);
                vp[t][i] = v[t][i] - m.alpha*(-dt*cr[t+1][i] + cv[t][i] - cv[t+1][i]);
                up[t][i] = u[t][i] - m.alpha*(u[t][i] - dt4*cr[t][i] - dt3*cr[t+1][i] - dt2*cv[t][i] - dt2*cv[t+1][i] + cu[t][i] - cu[t+1][i] + lu[t]*(i == 2));
                wp[t][i] = (t != tF)*(w[t][i] - m.alpha*(param.traj_wei_w*w[t][i] - cu1[t+1][i]));
            }
        }

        // Initial Conditions
        rp[0] = m.r0;
        vp[0] = {0,0,0};
        up[0] = projConeBall(up[0], param.syst_utheta, param.syst_umax);
        wp[0] = projBall(wp[0], param.syst_wmax);

        // Final Conditions
        for(long unsigned int t = m.tau+1; t < tF+1; t++) {
            rp[t] = m.rf;
            vp[t] = m.vf;
            up[t] = m.uf;
            wp[t] = m.wf;
        }

        // Projection of Primal Variable
        for(long unsigned int t = 1; t < m.tau+1; t++) {
            rp[t] = projCylinder(rp[t], m.cen, m.dir, m.len, m.rad);
            vp[t] = projBall(vp[t], param.syst_vmax);
            up[t] = projConeBall(up[t], param.syst_utheta, param.syst_umax);
            wp[t] = projBall(wp[t], param.syst_wmax);
        }
        
        // Compute PI feedback
        for(long unsigned int t = 0; t < r1.size(); t++) {
            for(long unsigned int i = 0; i < r1[t].size(); i++) {
                r1[t][i] = 2*r[t][i] - r[t][i];
                v1[t][i] = 2*v[t][i] - v[t][i];
                u1[t][i] = 2*u[t][i] - u[t][i];
                w1[t][i] = 2*w[t][i] - w[t][i];

                r[t][i] = (1-param.rho)*r[t][i] + param.rho*rp[t][i];
                v[t][i] = (1-param.rho)*v[t][i] + param.rho*vp[t][i];
                u[t][i] = (1-param.rho)*u[t][i] + param.rho*up[t][i];
                w[t][i] = (1-param.rho)*w[t][i] + param.rho*wp[t][i];
            }
        }

        // Update Dual Variables
        for(long unsigned int t = 0; t < tF; t++) {
            for(long unsigned int i = 0; i < cr[t].size(); i++) {
                crp[t+1][i] = cr[t+1][i] + m.beta*(r1[t+1][i] - r1[t][i] - dt*v1[t][i] - dt3*u1[t][i] - dt4*u1[t+1][i] - 0.5*pow(dt,2)*g[i]);
                cvp[t+1][i] = cv[t+1][i] + m.beta*(v1[t+1][i] - v1[t][i] - dt2*u1[t][i] - dt2*u1[t+1][i] - dt*g[i]);
                cup[t+1][i] = cu[t+1][i] + m.beta*(u1[t+1][i] - u1[t][i] - w1[t][i]);
            }
        }

        for(long unsigned int t = 0; t < cr.size(); t++) {
            for(long unsigned int i = 0; i < cr[t].size(); i++) {
                cr[t][i] = (1-param.rho)*cr[t][i] + param.rho*crp[t][i];
                cv[t][i] = (1-param.rho)*cv[t][i] + param.rho*cvp[t][i];
                cu[t][i] = (1-param.rho)*cu[t][i] + param.rho*cup[t][i];
            }
        }

        errOld = err;
        err = 0.0;
        for(long unsigned int t = 0; t < lu.size(); t++) {
            lup[t] = min(0.0, lu[t] + m.beta*(u1[t][2] - param.syst_umin));
            lu[t] = (1-param.rho)*lu[t] + param.rho*lup[t];
            err += pow(lu[t] - lu1[t], 2);
        }

        err += powSum(cr,cr1) + powSum(cv,cv1) + powSum(cu,cu1);
        err = min(errOld, err/(param.rho*m.beta*tF));

        int iter1 = iter - initIter;
        if(iter1 > 0) {
            if(iter1 == 1) {
                errOld = err;
            }
            else if(iter1 % freq == 0) {
                if(err < m.eps) {
                    ratio = 0;
                    break;
                }
                else {
                    ratio = err/errOld;

                    if(ratio > m.ratio) {
                        break;
                    }
                    errOld = err;
                }
            }
        }

    }

    return err;
}














inline vector<double> bisec(Parameters &param, vector<vector<double>> &pos, struct Trajectory &traj) {
    struct minTime m;
    const long unsigned int num = pos.size();
    vector<double> tauOPT(num);
    vector<double> tMin(num-1);
    vector<double> tMax(num-1);
    
    for(long unsigned int t = 0; t < num-1; t++) {
        tMin[t] = floor(2.0*param.corr_len[t]/(param.syst_dt*param.corr_vmax));
        tMax[t] = ceil(2.0*param.corr_len[t]/(param.syst_dt*param.corr_vmin));
    }
    
    vector<double> tauTest(num);
    for(int i = 0; i < num-1; i++) {
        tauTest[i+1] = tauTest[i] + tMax[i];
    }
    tauOPT = tauTest;

    const double tEPS = 2.5;
    m.ratio = 0.985;
    m.eps = 1E-3;
    double soltime = 0.0;
    for(long unsigned int n = 0; n < num-1; n++) {
        double ub = tMax[n];
        double lb = tMin[n];
        
        while(ub-lb > tEPS) {
            int delT = floor(0.5*(ub+lb));
            tauTest[n+1] = tauTest[n] + delT;
            for(int k = n+1; k < num-1; k++) {
                tauTest[k+1] = tauTest[k] + tMax[k];
            }

            Trajectory oldTraj = traj;
            chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
            double status = pipg_traj(param, m, traj, tauTest);
            chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

            soltime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
            if(status > 0) {
                ub = delT;  //if pipg converges, decrease ub
                tauOPT = tauTest;
            } 
            else {
                lb = delT;  //if pipg diverges, increase lb
                traj = oldTraj;
            }

            tauTest = tauOPT;
        }
    }
    traj.solTime = soltime;

    return tauOPT;
}


vector<double> runSim(struct Parameters &param, vector<vector<double>> pos) {
    //## Conversion to Corridors
    bool debug = false;
    const long unsigned int numCorridors = pos.size()-1;
    vector<vector<double>> cen(numCorridors, vector<double>(3,0));
    vector<vector<double>> dir(numCorridors, vector<double>(3,0));
    vector<double> len(numCorridors);
    for(long unsigned int i = 0; i < numCorridors; i++) {
        for(long unsigned int j = 0; j < cen[i].size(); j++) {
            cen[i][j] = 0.5*(pos[i+1][j] + pos[i][j]);
            dir[i][j] = 0.5*(pos[i+1][j] - pos[i][j]);
        }
        double n = norm(dir[i]);
        for(long unsigned int j = 0; j < cen[i].size(); j++) {
            dir[i][j] = dir[i][j]/n;
        }
        len[i] = n;
    }
    vector<double> wayptRad(numCorridors);
    wayptRad.assign(numCorridors, 0.2*0.81);
    
    param.corr_cen = cen;
    param.corr_dir = dir;
    param.corr_len = len;
    param.corr_rad = wayptRad;
    param.corr_vmax = 1.0*param.syst_vmax;    // maximum speed: m/s
    param.corr_vmin = 0.5*param.syst_vmax;    // minimum speed: m/s


    //## Final Constraints
    param.traj_r0 = pos.front();
    param.traj_rf = pos.back();


    //# TRIGGERING TIME CALCULATIONS    
    //## Quasiconvex Optimization
    Trajectory traj;
    chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    param.corr_time = bisec(param, pos, traj);
    chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    const double cost = (0.5*simplePow(traj.u) + 0.5*param.traj_wei_w*simplePow(traj.w))/(traj.u).size();

    if(debug) {
        cout << "\nnumCorridors = " << numCorridors << endl;
        cout << "PIPG Time = " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " [ms]" << endl;
        cout << "Cost = " << cost << endl << endl;
    }

    vector<double> statblock(4);
    statblock[0] = (chrono::duration_cast<chrono::milliseconds>(t2 - t1).count())/1000.0;
    statblock[1] = traj.solTime/1000.0; // soltime
    statblock[2] = -1; // parsertime
    statblock[3] = cost;

    return statblock;
}




int main() {
    struct Parameters param;
    param.syst_dt = 0.2;                  // sampling time    
    param.syst_m = 0.35;                  // mass of ACL quadrotor
    param.syst_g = {0.0,0.0,-9.81};       // gravity constant
    param.syst_umin = 2.0;                // minimum thrust magnitude
    param.syst_umax = 5.0;                // maximum thrust maginitude
    param.syst_vmax = 3.0;                // maximum speed
    param.syst_wmax = 3.0;                // maximum trust rate
    param.syst_utheta = M_PI/4.0;         // maximum tilting angle

    vector<double> uf_calc = param.syst_g;
    transform(uf_calc.begin(), uf_calc.end(), uf_calc.begin(),
                [param](double &g_i) {return -1.0*g_i*param.syst_m;});
    param.traj_uf = uf_calc;            // thrust for hovering
    param.traj_v0 = {0.0, 0.0, 0.0};    // inital velocity
    param.traj_vf = {0.0, 0.0, 0.0};    // final velocity
    param.traj_wf = {0.0, 0.0, 0.0};    // final thrust rate
    param.traj_wei_w = 1.0;             // weight on input rate (no larger than 1)

    param.maxIter = 10000;
    param.omega = 1.0;                  // ratio between primal and dual step sizes
    param.rho = 1.9;                    // extrapolation step size
    param.lam = max(1.0, param.traj_wei_w);

    
    
    // WAYPOINTS
    // Manual
    vector<vector<double>> pos = {{-4,-1,1}, {-2,1,1}, {0,-1,1}, {2,1,1}};
    
    // Variable
    // int numCorr = 3;
    // vector<vector<double>> pos(numCorr+1, vector<double>(3,0));
    // for(int n = 0; n < pos.size(); n++) {
    //     pos[n] = {-6.0+2.0*n, pow(-1.0,n), 1.0};
    // }

    // Conversion to Corridors
    const long unsigned int numCorridors = pos.size()-1;
    vector<vector<double>> cen(numCorridors, vector<double>(3,0));
    vector<vector<double>> dir(numCorridors, vector<double>(3,0));
    vector<double> len(numCorridors);
    for(long unsigned int i = 0; i < numCorridors; i++) {
        for(long unsigned int j = 0; j < cen[i].size(); j++) {
            cen[i][j] = 0.5*(pos[i+1][j] + pos[i][j]);
            dir[i][j] = 0.5*(pos[i+1][j] - pos[i][j]);
        }
        double n = norm(dir[i]);
        for(long unsigned int j = 0; j < cen[i].size(); j++) {
            dir[i][j] = dir[i][j]/n;
        }
        len[i] = n;
    }
    vector<double> wayptRad(numCorridors);
    wayptRad.assign(numCorridors, 0.2*0.81);
    
    param.corr_cen = cen;
    param.corr_dir = dir;
    param.corr_len = len;
    param.corr_rad = wayptRad;
    param.corr_vmax = 1.0*param.syst_vmax;    // maximum speed: m/s
    param.corr_vmin = 0.5*param.syst_vmax;    // minimum speed: m/s


    //## Final Constraints
    param.traj_r0 = pos.front();
    param.traj_rf = pos.back();


    //# TRIGGERING TIME CALCULATIONS    
    //## Quasiconvex Optimization
    Trajectory traj;
    chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    param.corr_time = bisec(param, pos, traj);
    chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    const double cost = (0.5*simplePow(traj.u) + 0.5*param.traj_wei_w*simplePow(traj.w))/(traj.u).size();

    cout << "\nnumCorridors = " << numCorridors << endl;
    cout << "PIPG Time = " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " [ms]" << endl;
    cout << "Cost = " << cost << endl << endl;

    return 0;
}