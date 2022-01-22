# Acoustic-SCR
Acoustic SCR MCMC samplers

Currently, only the case of stationary callers and known individual ID is in this repository. This is the same model as described by Stevenson et al (2021), but fit via MCMC. The current assumption is that the call rate is Poisson, but it can be changed to any count distribution (but zero truncation requires modifying the custom call number update to prevent proposing that an individual may have 0 calls).

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13522

See testscript.