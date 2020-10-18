# Stan Selection Models
 
This repository contains an initial attempt to implement Bayesian selection models in Stan. This repository is superseded by the RoBMA package (https://github.com/FBartos/RoBMA) that overcomes many of the problems.

### The main issues with the Stan implementation:
- the HMC does not deal well with the step likelihood function due to the discontinuities in the likelihood.
- the sampling of non-step functions is extremely slow due to the need to numerically integrate the denominator of the likelihood function.

#### The implementation was not properly tested and might contain errors.