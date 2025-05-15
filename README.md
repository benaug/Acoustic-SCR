# Acoustic-SCR
Acoustic SCR MCMC samplers

This repository contains MCMC samplers for acoustic SCR that I am converting over to nimble. I consider stationary and mobile callers (bivariate normal movement) and known, unknown, and partial individual IDs. I have provided examples here for stationary and mobile callers with known individual IDs. I will upload the remainder after a preprint is posted. The slides from my Euring 2021 talk are in this repository. Those results are from the R-based samplers, nimble should be faster. Also, I truncated the BVN cluster dispersion model by the state space extent there, but I don't do that in nimble because it is slower and not necessary. You just need to make sure you buffer the ARU array by at least 1/2 the 95-99% BVN radius or so.

The model for stationary callers with known ID is the same as the model described by Stevenson et al (2021), but fit via MCMC and with an attenuation function observation model. For this model, there are 2 testscripts, one assuming that the call rate is Poisson and another assuming it is zero-truncated negative binomial. Any count distribution is allowed, but care should be taken that the custom count update is correct for a given count model.

The MCMC sampler for mobile callers is more complicated and is very "hard coded", meaning that many modifications to the BUGS code will require changes to the custom Metropolis-Hastings updates. Also, it requires some tuning by the user. Use and especially modify at your own risk. A novel aspect of this approach is that it uses double data augmentation, allowing for individual covariates and call covariates. This approach may be generally useful for cluster SCR problems where both "parents" and "children" each have individual covariates.

Unsolicited opinion: I think that a zero-truncated count model is most appropriate for acoustic SCR and we should estimate the "calling population size" because the population size of 0 callers is not identifiable if it differs from that expected by a count distribution that is not zero-truncated. E.g., if the call distribution is Poisson with lambda=10, we expect essentially no individuals that do not call during the survey (dpois(0,10)), though the component of the population that does not call could be vastly different. If so, the 0 caller component of the population is not estimable without more assumptions and replication (e.g., multiple occasions with an individual-level call/no call process). An alternative approach is that you could use a count distribution whose support includes 0, and derive the calling population size by summing the individuals in the model that currently have 1 or more counts.

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13522

The attenuation function is from Efford and Dawson (2009) and Dawson and Efford (2009)

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/08-1735.1

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2009.01731.x

These models use count prior data augmentation for call number only. Could be used for number of individuals as well. 
https://github.com/benaug/SCR-Count-Prior-Data-Augmentation