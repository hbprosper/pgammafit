# pgammafit
Implementation of Bayesian fit using the Poisson-gamma model, based on the paper at https://inspirehep.net/record/441459

### Notes

* This package depends on Root  and the root extended mathematical
library MathMore, which, in turn, depends on the GNU Scientific
Library (GSL). Therefore, GSL must be available on your system so that
configure script of Root can find it, otherwise your build of Root
will not contain MathMore.

* Example
```
cd example
python pgammafit.py
```
