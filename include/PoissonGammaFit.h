#ifndef POISSONGAMMAFIT_H
#define POISSONGAMMAFIT_H
///////////////////////////////////////////////////////////////////////////////
// File:    PoissonGammaFit.h
// Description:	
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
//              f is a vector of vectors of size N x M scale factors.
//              such that f*A = Aeff, f*dA = sqrt(Aeff), where Aeff
//              is an effective Poisson count.
//
// Created: 27-Jun-2005 Harrison B. Prosper
//          19-Feb-2011 HBP add add(..)
//          16-Apr-2011 HBP allow for parameter fixing
//          05-Jun-2016 HBP allow different scaling assuptions for different
//                      sources.
///////////////////////////////////////////////////////////////////////////////
#include <vector>

typedef std::vector<bool> vbool;
typedef std::vector<double>  vdouble;
typedef std::vector<std::vector<double> > vvdouble;

namespace pg {
  double 
  poissongamma(vdouble&	      D,   // Counts  "D_i" for data.
	       vdouble&	      p,   // Weights "p_j" 
	       vvdouble&      A,   // Counts  "A_ji" for up to 8 sources
	       vvdouble&      f,   // scale factor for  "A_ji"
	       vbool& scale);      // Scale p_j if true
}

/** Fit Poisson/gamma model to data.
    The data consists of a sequence of \f$i = 1 \ldots M\f$ counts 
    \f$D_i\f$ (i.e., a histogram),
    presumed to arise from a mixture of sources, say \f$t\bar{t}\f$,
    \f$W + jets\f$, and \f$QCD\f$ events. The goal is to find the yield of
    each source, given \f$j = 1 \ldots N\f$ sequences of counts 
    \f$A_{ji}\f$ that describe the ``shapes" of the \f$N\f$ sources.
    In other words, given a data histogram \f$D\f$ and \f$N\f$ histograms 
    \f$A_j = A_{j1}, \ldots, A_{jM}\f$, one for each source \f$j\f$, 
    we wish to find the yield of each. This class solves the
    problem using a Bayesian approach \ref bayes1 "[1]".
<p>
    <b>Model</b>
    <p>
    - <i>Likelihood</i> - The likelihood function is assumed to be given by
    \f[ 
    p(D|p, a) = \prod_{i=1}^M\prod_{j=1}^N \mbox{Poisson}(D_i|d_i),
    \f]
    where
    <br>
    \f$ d_i = \sum_{j=1}^N \, p_j \,  a_{ji},  \, \, \, i=1 \ldots M.\f$
    <br>
    The unknown parameters \f$p_j\f$ and \f$a_{ji}\f$ are 
    the source ``strengths" and
    expected source counts, respectively.
    <p>
    - <i>Prior</i> - The prior is factorized as follows 
    \f$
    \pi(p, a) = \pi(a|p) \, \pi(p)
    \f$ and 
    \f$\pi(a|p)\f$ is assumed to be a product of gamma densities:
    \f[ p(a|p) = \prod_{i=1}^M\prod_{j=1}^N 
    \mbox{gamma}(a_{ji}|A_{ji}+1, 1),
    \f]
    where
    \f[
    \mbox{gamma}(x|\alpha,\beta) = 
    \beta^{\alpha} x^{\alpha-1} e^{-\beta x} / \Gamma(\alpha).
    \f]
    <p>
    - <i>Posterior</i> - The posterior density of the source ``strengths" is
    given by 
    \f[
    p(p|D) = p(D|p) \, \pi(p) \, / \, \int p(D|p) \, \pi(p) \, dp,
    \f]
    where the marginal likelihood \f$p(D|p)\f$, that is, \f$p(D|p, a)\f$
    integrated with respect to the nuisance parameters \f$a_{ji}\f$,
    is given in Ref.\ref bayes1 "[1]".
    
    The form of the prior follows from the assumption that each count 
    \f$A_{ji}\f$ is associated with a Poisson distribution and that the
    the magnitude of the count \f$A_{ji}\f$ is such that the 
    likelihood dominates
    the prior. When that assumption is not justified, a reasonable 
    procedure is to use
    Jeffreys prior\ref bayes2 "[2]",
which leads to the replacement of
    \f$A_{ji}\f$ with \f$A_{ji}-1/2\f$. For most practical applications, 
    however, this will not make much of a difference. 
    <p>
    <b>Weighted histograms</b>
    <p>
    For weighted histograms, the Baysian fit requires the specification of
    the uncertainty \f$\delta A_{ji}\f$ associatd with the estimate 
    \f$A_{ji}\f$ since we can no longer assume the uncertainty is 
    \f$\sqrt{a_{ji}}\f$. 
    Internally, each count \f$A_{ji}\f$ is scaled by the  
    factor
    \f$f_{ji} = (A_{ji} / \delta A_{ji} )^2\f$ so that the 
    relative uncertainty in
    the effective bin count \f$A_{ji, eff} = f_{ji} A_{ji}\f$ matches 
    \f$\delta A_{ji} / A_{ji}\f$. In the limit 
    \f$f_{ji} \rightarrow \infty\f$, as 
    expected, the marginal likelihood becomes the same as the likelihood with
    \f$a_{ji}\f$ replaced by \f$A_{ji}\f$.

<p>
<b>References</b>
<br>
\anchor bayes1 [1] P.C. Bhat, H.B. Prosper, and 
S.S. Snyder, <i>Bayesian analysis of multisource data</i>, 
Phys. Lett. B<b>407</b>:73-78, 1997, 
<a href="http://lss.fnal.gov/archive/test-preprint/fermilab-pub-96-397.shtml">
FERMILAB-PUB-96-397</a>.
<br>
\anchor bayes2 [2] L. Demortier, S. Jain, and H.B. Prosper,
    <i>Reference priors for high energy physics</i>, 
    Phys. Rev. D<b>82</b>:034002, 2010, 
<a href="http://arxiv.org/pdf/1002.1111">e-Print: arXiv:1002.1111</a>.
 */
class PoissonGammaFit
{
 public: 

  /// Status codes for fit
  enum Status
    {
      kSUCCESS  = 0,
      kNMISMATCH=-1,
      kMMISMATCH=-2,
      kZEROBINS =-3,
      kZEROSRCS =-4,
      kTOOMANYSRCS=-5,
      kNOTCONVERGED=-6
    };

  ///
  PoissonGammaFit();

  ///
  PoissonGammaFit(vdouble&	D,  // Data vector
                  int  verbosity=0);
  ///
  virtual ~PoissonGammaFit();

  /// True if all is well.
  bool       good();

  /** Add a source with a floating yield by default.
      @param A  - source counts
      @param dA - source uncertainties
      @param scale - if true, scale p for this source so that p is a count.
      @param fixed - if true, keep this parameter fixed.
   */
  void       add(vdouble& A, vdouble& dA, bool scale=true, bool fixed=false);

  /// Execute the fit.
  bool       execute(std::vector<double>& guess);

  /// Return status code.
  Status     status();

  /// Return mode of posterior density
  vdouble&   mode();

  /// Return mean of posterior density
  vdouble&   mean();
    
  /// Return width of posterior density. 
  vdouble&   width();

  /// Return initial guesses of parameters.
  vdouble&   guess();

  /// Return which parameters are fixed, if any
  std::vector<bool>&   fixed();

  /// Return error matrix of parameters
  vvdouble&  errorMatrix();

  /// Return logarithm of evidence.
  double     logEvidence();

  /// Return logarithm of likelihood at its mode.
  double     logLikelihoodMax();

  ///
  double     logLikelihood(vdouble& x);

  ///
  double     logLikelihood(double* x);

 private:

  vdouble  _D;
  vvdouble _A;
  vvdouble _f;
  vvdouble _a;
  vbool    _scale;  
  vbool    _fixed;
  vdouble  _guess;
  vdouble  _ns;
  vdouble  _mode;
  vdouble  _mean;
  vdouble  _width;
  vvdouble _cov;
  
  int      _verbosity;  
  Status   _status;
  int      _N;
  int      _M;
  double   _nd;
  double   _loglikemax;
  double   _logevidence;

  bool     _findmode();

/* #ifdef __WITH_CINT__ */
/*  public: */
/*   ClassDef(PoissonGammaFit, 1) */
/* #endif */
};


#endif
