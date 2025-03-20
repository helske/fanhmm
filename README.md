Codes for reproducing the results of the paper *Feedback-augmented Non-homogeneous Hidden Markov Models for Longitudinal Causal Inference*
================

- Folder `application` contains codes and results relating to the parental leave reform application, as well as the synthetic data `synthetic_leave_data.rds` and corresponding model fit.
- Folder `figures` contains all figures of the paper and the supplement, as well as two additional figures `ace_S3_reform_states.pdf` and `ace_S3_reform_conditionals.pdf` which show the causal effect of the reform on state marginals and conditional emissions respectively.
- Folder `simulation_experiments` contains all codes and simulations used in the simulation experiments of the main paper and the supplement.

All the modeling codes require the `seqHMM` R package, which can be installed from the source package in the repo using `install.packages("seqHMM_2.0.0.tar.gz")` after downloading the package. The latest version of `seqHMM` is available at https://github.com/helske/seqHMM/.
