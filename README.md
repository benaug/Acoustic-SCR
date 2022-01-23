# Acoustic-SCR
Acoustic SCR MCMC samplers

Currently, only the case of stationary callers and known individual ID is in this repository. This is the same model as described by Stevenson et al (2021), but fit via MCMC and with an attenuation function observation model. 

There are currently 2 testscripts, one assuming that the call rate is Poisson and another assuming it is zero-truncated negative binomial. Any count distribution is allowed, but care should be taken that the custom count update is correct for a given count model.

My opinion is that a zero-truncated count model is most appropriate and we should estimate the "calling population size" because the population size of 0 callers is not identifiable if it differs from that expected by a count distribution that is not zero-truncated. E.g., if the call distribution is Poisson with lambda=10, we expect essentially no individuals that do not call during the survey (dpois(0,10)), though the component of the population that does not call could be vastly different. If so, the 0 caller component of the population is not estimable without more assumptions and replication (e.g., multiple occasions with an individual-level call/no call process).

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13522

The attenuation function is from Efford and Dawson (2009) and Dawson and Efford (2009)

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/08-1735.1

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2009.01731.x

See testscript.