# SW-CCR
Balance metric for use in covariate-constrained randomization in stepped-wedge cluster randomized trials.

A method for covariate-constrained randomization in cluster randomized trials (CRT) or stepped wedge (SW) designs involves assessing
the balance of covariates for all possible randomizations. The randomization scheme used for the trial is selected from a restricted
pool of possible randomizations with acceptable balance based on a defined metric. Our proposed metric ("B") for SW designs is a modified
version of a metric commonly used in CRT. The attached code calculates "B" for all possible randomizations of a user-specified SW design
and covariate z-scores. It selects one randomization scheme from a pool of those with the best balance (user-specified percentile cutoff 
of the distribution of B).
