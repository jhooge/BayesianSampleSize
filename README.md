# Bayesian Sample Size Calculation

This repository contains code for a shiny application to estimate sample sizes 
and visualize the model for clinical trials.
It assumes a beta-binomaial model and a single treatment class. Based on the number of successes 
(drawn from a binomial distribution) and prior knowledge (following a beta distribution),
the number of samples can be computed.

![A fancy .gif presenting the app](img/BayesianSampleSize.gif)

## Motivation
In clinical research, parameters required for sample size calculation are 
usually unknown. A typical approach is to use estimates from some pilot 
studies as the true parameters in the calculation. 
This approach, however, does not take into consideration sampling error. 
Thus, the resulting sample size could be misleading if the sampling error 
is substantial. As an alternative, we suggest a Bayesian approach to include
prior knowledge about the underlying 

## Author
Jens Hooge [jens.hooge@bayer.com](mailto:jens.hooge@bayer.com)

