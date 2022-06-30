#include <vector>
using namespace std;


// STRUCTS
struct Parameters;
struct Trajectory;
struct minTime;

// HELPER FUNCTIONS

// Inputs = c (input vector of size n)
// Outputs = a (double)
double norm(vector<double> &c);

// Inputs = u (input vector of size n)
// Outputs = printline
void printVector(vector<double> &u);

// Inputs = A (input matrix of size NxM)
// Outputs = printline
void printMatrix(vector<vector<double>> &A);

// Inputs = A (input matrix of size NxM), x (input vector of size n)
// Outputs = v (input vector of size n)
vector<double> matVec(vector<vector<double>> &A, vector<double> &x);

// Inputs = A (input matrix of size NxM)
// Outputs = T (input matrix of size MxN)
vector<vector<double>> transpose(vector<vector<double>> &A);

// Helper Function for Algo 6: Norm Power
// Inputs = x (input matrix of size NxM),
//          u (input matrix of size NxM)
// Outputs = s (output double)
double normPower(vector<vector<double>> &x, vector<vector<double>> &u);

// Helper Function for Algo 7: normDiff
//// For now, assume dimensions correct (.size=row, i.size=col)
// Inputs = v (input matrix of size),
//          vBar (input matrix of size),
// Outputs = s (output double)
double normDiff(vector<vector<double>> &v, vector<vector<double>> &vBar);


double powSum(vector<vector<double>> &x, vector<vector<double>> &x1);


double simplePow(vector<vector<double>> &x);



// PROJECTION FUNCTIONS

// Algo 1: Projection onto the nonnegative orthant
// Inputs = c (input vector of size n)
// Outputs = x (output vector of size n)
vector<double> projPosOrth(vector<double> &c);

// Algo 2: Projection onto a box
// Inputs = c (input vector of size n), 
//          l (lower bound as a real number), 
//          u (upper bound as a real number)
// Outputs = x (output vector of size n)
vector<double> projBox(vector<double> &c, vector<double> &l, vector<double> &u);

// Algo 3: Projection onto a ball
// Inputs = c (input vector of size n), 
//          p (ball radius as a positive real number, usually rho)
// Outputs = x (output vector of size n)
vector<double> projBall(vector<double> &c, double &p);

// Algo 4: Projection onto an (icecream) cone
// Inputs = c (input vector of size n), 
//          g (cone ratio as a positive real number, usually gamma)
// Outputs = x (output vector of size n)
vector<double> projCone(vector<double> &c, double &g);

// Algo 5: Projection onto the intersection of an icecream cone and a ball
// Inputs = c (input vector of size n), 
//          g (cone ratio as a positive real number, usually gamma), 
//          p (ball radius as a positive real number, usually rho)
// Outputs = x (output vector of size n)
vector<double> projConeBall(vector<double> &c, double &g, double &p);

// Algo 6: Projection onto a cylinder
// Inputs = c (input vector of size n), 
//          eta (a positive real number, usually eta), 
//          p (ball radius as a positive real number, usually rho), 
//          h (vector of size n), 
//          e (vector of size n with norm = 1)
// Outputs = x (output vector of size n)
vector<double> projCylinder(vector<double> &c, vector<double> &h, vector<double> &e, double &eta, double &p);




// LARGER FORM ALGORITHMS

// new power_iter function translation to C++
// Helper function for pipg_traj
double power_iter(const long unsigned int &tau, struct Parameters &param);


// pipg_traj C++ translation
// Inputs = param (struct)
// Outputs = 
struct Trajectory pipg_traj(struct Parameters &param);



double pt2pt(struct Parameters &param, struct minTime &m);



vector<double> bisec(Parameters &param, vector<vector<double>> &pos);