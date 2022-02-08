# RandomThinIDCov

Nimble samplers for the random thinning SCR model from Jimenez et al. (2021).

https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7091

This repo houses samplers that combine the random thinning model with categorical SCR from Augustine et al. (2019).

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2627

As a result, 

1) The individual IDs of unidentified samples are probabilistically resolved with less uncertainty
2) Model parameters can be a function of the (partial) individual ID covariates. This includes typical SCR parameters like lam0 and sigma, but also the thinning rates. E.g., male cats may have a larger sigma than female cats and calico cats are easier to identify to individual than black cats (substitute any categorical covariates here.)


(only base RT model uploaded so far)