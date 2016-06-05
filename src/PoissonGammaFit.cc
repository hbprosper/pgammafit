//-----------------------------------------------------------------------------
// File:	PoissonGammaFit.cc
// Description:	Fit up to 8 sources to data.
//
//              Model:
//
//              d_i = Sum p_j * a_ji,  i=1..M,   j=1..N
// 
//              D_i  is count corresponding to true mean d_i
//              A_ji is count corresponding to true mean a_ji
//
//              D is a vector of M observed counts (that is, over M bins)
//              p is a vector of N parameters (that is, over N sources)
//              A is a vector of vectors of size N x M counts.
//
//              Method:
//
//              The integration over source priors is done exactly, so beware
//              of histograms with huge counts!
//
//              The logarithm of the evidence is approximated by the expression
//              (see for example, Probability Theory: The Logic of Science,
//              page 2404, by the late Prof. Edwin Jaynes)
//
//              log(p(D|M,I)) ~ log(Lmax) + (N/2)*log(2*pi) 
//                            + (1/2) log(det(Sigma))
//
// Created:     27-Jun-2005 Harrison B. Prosper
// Updated:     06-Feb-2006 HBP Adapt for use in Root
//              19-Feb-2011 HBP generalize to handle weighted events
//              16-Apr-2011 HBP allow some parameters to be fixed
//              03-Jun-2016 HBP migrate to github
//-----------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "PoissonGammaFit.h"
#include "TMinuit.h"
#include "TMatrix.h"

// #ifdef __WITH_CINT__
//   ClassImp(PoissonGammaFit)
// #endif

using namespace std;
//-----------------------------------------------------------------------------

static PoissonGammaFit* fitobject=0;
namespace
{
  void logLike(int&    /** npar */, 
               double* /** grad */, 
               double& fval,
               double* xval,
               int    /**iflag */)
  {
    vector<bool>& fixed = fitobject->fixed();
    vector<double>& guess = fitobject->guess();
    vector<double> p(fixed.size());
    int j=0;
    for(unsigned int i=0; i < fixed.size(); ++i)
      {
        if ( fixed[i] )
          {
            p[i] = guess[i]; // use user supplied value
          }
        else
          {
            p[i] = xval[j];
            j++;
          }
      }
    fval  = -fitobject->logLikelihood(p);
  }

  const int MAXITER=10000;
  const int MAXSRCS=8; // Maximum number of sources
}

//-----------------------------------------------------------------------------
// Class to perform Bayesian fit
//-----------------------------------------------------------------------------
PoissonGammaFit::PoissonGammaFit()
  : _D(vdouble()),
    _A(vvdouble()),
    _f(vvdouble()),
    _a(vvdouble()),
    _scale(vbool()),
    _fixed(vbool()),
    _guess(vdouble()),
    _ns(vdouble()),
    _mode(vdouble()),
    _mean(vdouble()),
    _width(vdouble()),
    _cov(vvdouble()),
    _verbosity(0),
    _status(kSUCCESS),
    _N(0),
    _M(0),
    _nd(0),
    _loglikemax(0),
    _logevidence(0)
{}

PoissonGammaFit::PoissonGammaFit(vdouble& D,  // Counts  "D_i" for data.
                                 int verbosity)
  : _D(D),
    _A(vvdouble()),
    _f(vvdouble()),
    _a(vvdouble()),
    _scale(vbool()),
    _fixed(vbool()),
    _guess(vdouble()),
    _ns(vdouble()),
    _mode(vdouble()),
    _mean(vdouble()),
    _width(vdouble()),
    _cov(vvdouble()),
    _verbosity(verbosity),
    _status(kSUCCESS),
    _N(0),
    _M(0),
    _nd(0),
    _loglikemax(0),
    _logevidence(0)
{
  if ( _verbosity > 0 )
    cout << "PoissonGammaFit - Started" << endl;

  _status = kSUCCESS;
 
 _M = _D.size();  // Number of bins    (M)

  if ( _M < 1 )
    {
      std::cout 
	<< "**Error - PoissonGammaFit - bin count is zero"
	<< std::endl;
      _status = kZEROBINS;
      return;      
    }
 
  // Get data count

  _nd = 0.0; 					
  for ( int i = 0; i < _M; i++ ) _nd += _D[i];
}

void 
PoissonGammaFit::add(vector<double>& A, vector<double>& dA,
		     bool scale,
		     bool isfixed)
{
  // Compute scale factors
  vector<double> f(A.size(), 1); // default scale factor = 1
  if ( dA.size() == A.size() )
    {
      for(unsigned int i=0; i < A.size(); ++i)
        {
          if ( dA[i] > 0 )
            {
	      // assume:
	      // f *  A = Aeff
	      // f * dA = sqrt(Aeff)
	      // so f = A / (dA*dA)
	      // where Aeff is the effective count
              f[i] = A[i] / (dA[i]*dA[i]);
            }
        }
    }

  _A.push_back(A);
  _a.push_back(A);
  _scale.push_back(scale);
  _fixed.push_back(isfixed);
  _f.push_back(f);
}

bool
PoissonGammaFit::execute(vector<double>& someguess)
{
  _N = _A.size();  // Number of sources (N)  
  if ( _N < 1 )
    {
      std::cout 
        << "**Error - PoissonGammaFit - zero sources"
        << std::endl;
      _status = kZEROSRCS;
      return false;
    }

  if ( _N > MAXSRCS )
    {
      std::cout 
        << "**Error - PoissonGammaFit - too many sources: " << _N 
        << std::endl;
      _status = kTOOMANYSRCS;
      return false;
    }

  if ( _A[0].size() != (unsigned int)_M )
    {
      std::cout 
        << "**Error - PoissonGammaFit - mis-match in number of sources\n"
        << "size(D): " << _M << " differs from size(A[0]) = " 
        << _A[0].size()
        << std::endl;
      _status = kMMISMATCH;
      return false;
    }

  // Get total counts per source
  _cov.clear();
  _mode.clear();
  _mean.clear();
  _width.clear();
  _ns.clear();
  _guess.clear();

  for ( int j = 0; j < _N; j++ )
    {
      _cov.push_back(vector<double>(_N));
      _mode.push_back(0);
      _mean.push_back(0);
      _width.push_back(0);
      _ns.push_back(0);
      _guess.push_back(someguess[j]);

      for ( int i = 0; i < _M; i++ ) _ns[j] += _A[j][i];
    }

  // Ok, use Minuit to find the mode of marginal likelihood
  return _findmode();
}

PoissonGammaFit::~PoissonGammaFit() {}

double 
PoissonGammaFit::logLikelihood(vdouble& p)
{
  for(int i=0; i < _N; i++) if (p[i] < 0) return 0;
  return pg::poissongamma(_D, p, _A, _f, _scale);
}

double 
PoissonGammaFit::logLikelihood(double* point)
{
  vdouble p(_N);
  copy(point, point+_N, p.begin());
  return pg::poissongamma(_D, p, _A, _f, _scale);
}

bool 
PoissonGammaFit::good() { return _status == kSUCCESS; }

PoissonGammaFit::Status
PoissonGammaFit::status() { return _status; }

vdouble&
PoissonGammaFit::mode() { return _mode; }

vdouble&
PoissonGammaFit::mean() { return _mean; }

vdouble&
PoissonGammaFit::width() { return _width; }

vdouble&
PoissonGammaFit::guess() { return _guess; }

vector<bool>&
PoissonGammaFit::fixed() { return _fixed; }

vvdouble&
PoissonGammaFit::errorMatrix() { return _cov; }

double 
PoissonGammaFit::logEvidence() { return _logevidence; }

double 
PoissonGammaFit::logLikelihoodMax() { return _loglikemax; }

bool
PoissonGammaFit::_findmode()
{
  fitobject = this; // Update global pointer to this object

  _status = kSUCCESS;

  if ( _verbosity > 0 )
    cout << "\tFind mode" << endl;

 // Count number of free parameters
  int N = 0;
  for(int i=0; i < _N; ++i) if ( !_fixed[i] ) N++;

  // Find mode of likelihood function

  TMinuit minuit(N);
  minuit.SetFCN(logLike);

  int ierflag;
  int j=0;
  for (int i=0; i < _N; i++)
    {
      if ( !_fixed[i] )
        {
          char name[80]; sprintf(name,"x%d", j);
          double x    = _guess[i];
          double step = x / 1000;
          double minx = 0;
          double maxx = 10 * x;
          minuit.mnparm(j, name, x, step, minx, maxx, ierflag);
          j++;
        }
    }

  // Do fit

  minuit.SetErrorDef(0.5);
  minuit.SetMaxIterations(MAXITER);
  minuit.Migrad();

  // Check for convergence

  double fmin, fedm, errdef;
  int nvpar, nparx, istatus;
  minuit.mnstat(fmin, fedm, errdef, nvpar, nparx, istatus);
  if ( istatus == 3 )
    {
      // Get parameters
      _loglikemax = -fmin;
      int k=0;
      for (int i=0; i < _N; i++)
        {
          if ( _fixed[i] )
            {
              _mode[i] = _guess[i];
              _width[i]= 0;
            }
          else
            {
              minuit.GetParameter(k, _mode[i], _width[i]);
              k++;
            }
        }
    }
  else
    {
      if ( _verbosity > 0 )
	  cout << "\tMode finding FAILED" << endl;
      _status = kNOTCONVERGED;
      return false;
    }

  if ( _verbosity > 0 )
    cout << "\tMode FOUND" << endl;
  
  // Compute estimate of log(evidence)
  // log p(D|M,I) ~ ln Lmax + (N/2) ln 2*pi + (1/2) ln det(Sigma)
  
  int ndim = N * N;
  double errmat[2500];
  minuit.mnemat(&errmat[0], ndim);

  TMatrix mat(N, N);
  for(int ii=0; ii < N; ++ii)
    for(int jj=0; jj < N; ++jj)
      { 
        _cov[ii][jj] = errmat[ii*ndim+jj];
        mat[ii][jj] = _cov[ii][jj];
      }
  _logevidence = _loglikemax + 0.5*N * log(2*M_PI);
  double det = mat.Determinant();
  if ( det > 0 ) _logevidence += log(det);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
// File:	poissongamma.cc
// Description:	x = poissongamma(D,p,A,f)
// 
//              compute the marginal likelihood:
//              a product of Poisson distributions
//              and a prior that is a product of gamma distributions.
//
//              Model:
//
//              d_i = Sum p_j * a_ji,  i=1..M,   j=1..N
// 
//              D_i  is count corresponding to true mean d_i
//              A_ji is count corresponding to true mean a_ji
//              f_ji is an optional scale factor associated with A_ji
//
//              D is a vector of M observed counts (that is, over M bins)
//              p is a vector of N parameters (that is, over N sources)
//              A is a vector of vectors of size N x M counts.
//
//              Simple 1-D Application:
//              
//              d =  xsec * e_a * a + e_b * b 
//
//              where p_1 = xsec * e_a
//                    p_2 = e_b
//                    a_11= a
//                    a_21= b
//
//              and e_a, e_b are scale factors defined by
//              e = sigma^2 / estimate
//
// WARNING: The number of terms grows exponentially with the number of 
//          sources and the maximum bin count. Therefore, this routine
//          really only makes sense for rather few sources with modest
//          bin counts. Even so, we allow counts per bin up to 50,000.
//
// Created:     20-Dec-2001 Harrison B. Prosper
//              11-Mar-2004 HBP add additional interfaces
//              07-Jun-2005 HBP increase number of sources to 10! See warning
//                          above.
//              19-Feb-2011 HBP generalize to make use of user-supplied
//                          scale factors. This is needed for weighted 
//                          histograms. Reduce maximum number of sources to 6
//              27-Apr-2011 HBP increase number of sources to 8
///////////////////////////////////////////////////////////////////////////////


// Global variables to avoid memory allocation.
namespace pg {

  const int MAXSRC = MAXSRCS;     // Maximum number of sources
  const int MAXBUF = 50000;       // Maximum count per bin

  static long double c[MAXSRC][MAXBUF];
  static double ns[MAXSRC];
  static double s [MAXSRC];
  static double x [MAXSRC];
  static double y [MAXSRC];

  double 
    poissongamma(vdouble&	D,   // Counts  "D_i" for data.
                 vdouble&	p,   // Weights "p_j" 
                 vvdouble&	A,   // Counts  "A_ji" for up to 8 sources
                 vvdouble&	f,   // scale factor for  "A_ji"
                 vector<bool>& scale)     // Scale p_j if true  
  {
    int N = p.size(); // Number of sources (N)
    int M = D.size(); // Number of bins    (M)

    // Check inputs
    if ( A.size() != (unsigned int)N )
      {
        std::cout << "**Error - poissongamma - mis-match in number of sources"
                  << endl
                  << "size(p): " << N << " differs from size(A) = " << A.size()
                  << std::endl;
        exit(0);
      }

    if ( A[0].size() != (unsigned int)M )
      {
        std::cout << "**Error - poissongamma - mis-match in number of binss\n"
                  << "size(D): " << M << " differs from size(A[0]) = " 
                  << A[0].size()
                  << std::endl;
        exit(0);
      }

    if ( M < 1 ) return -1.0;
  
    if ( ( N < 2 ) || ( N > MAXSRC ) ) return -2.0;
  
    // Get total counts per source

    for (int j = 0; j < N; ++j)
      {
        ns[j] = 0;
        for (int i = 0; i < M; ++i)
          ns[j] += A[j][i];
      }

    // loop over the M terms of the product,
    // corresponding to the M bins

    double prob=0;
    for (int i = 0; i < M; ++i)
      {
        int Di = (int)D[i]; // data count for bin i

        // compute terms of sum from zero to D

        // first do zero...      
        for (int j = 0; j < N; ++j)
          {
            // Optionally, normalize this source to unit area so that
            // x[i] becomes the actual source count.
            if ( scale[j] )
              x[j] = p[j] / (ns[j]+M);
            else
              x[j] = p[j];
            s[j] = A[j][i];

            // Apply user supplied scale factor
            // This is needed to take account of weighted histograms
            if ( (int)f[j].size() == M )
              {
                if ( f[j][i] > 0 )
                  {
                    x[j] /= f[j][i];
                    s[j] *= f[j][i];  
                  }
              }

            y[j] = x[j] / (1 + x[j]);
            c[j][0] = pow(1 + x[j], -(s[j] + 1));
          }
      
        // ...then 1 to D
        if ( Di > 0 )
          {
            for (int k = 1; k < Di+1; ++k)
              for (int j = 0; j < N; ++j)
                c[j][k] = c[j][k-1] * y[j] * (s[j] + k)/k;
          }

        // compute sum

        double sum = 0.0;
        switch (N)
          {
          case 2:
            for (int j = 0; j < Di+1; ++j)
              sum += 
                c[0][j] *
                c[1][Di-j];
            break;

          case 3:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                sum += 
                  c[0][j] *
                  c[1][k] *
                  c[2][Di-j-k];
            break;

          case 4:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  sum += 
                    c[0][j] *
                    c[1][k] *
                    c[2][l] *
                    c[3][Di-j-k-l];
            break;
          
          case 5:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    sum += 
                      c[0][j] *
                      c[1][k] *
                      c[2][l] *
                      c[3][m] *
                      c[4][Di-j-k-l-m];
            break;
          
          case 6:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      sum += 
                        c[0][j] *
                        c[1][k] *
                        c[2][l] *
                        c[3][m] *
                        c[4][n] *
                        c[5][Di-j-k-l-m-n];
            break;
          
          case 7:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                        sum += 
                          c[0][j] *
                          c[1][k] *
                          c[2][l] *
                          c[3][m] *
                          c[4][n] *
                          c[5][jj] *
                          c[6][Di-j-k-l-m-n-jj];
            break;

          case 8:
            for (int j = 0; j < Di+1; ++j)
              for (int k = 0; k < Di+1-j; ++k)
                for (int l = 0; l < Di+1-j-k; ++l)
                  for (int m = 0; m < Di+1-j-k-l; ++m)
                    for (int n = 0; n < Di+1-j-k-l-m; ++n)
                      for (int jj = 0; jj < Di+1-j-k-l-m-n; ++jj)
                        for (int kk = 0; kk < Di+1-j-k-l-m-n-jj; ++kk)
                          sum += 
                            c[0][j] *
                            c[1][k] *
                            c[2][l] *
                            c[3][m] *
                            c[4][n] *
                            c[5][jj] *
                            c[6][kk] *
                            c[7][Di-j-k-l-m-n-jj-kk];
            break;
          };
	prob += log(sum);
      }
    return prob;
  }
};
