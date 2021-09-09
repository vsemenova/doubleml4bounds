### Double Machine Learning for Set-Identified Linear Models

The code in this replication package replicates the simulation results contained in Tables 1,2 for the paper

Semenova, Vira "Debiased Machine Learning of Set-Identified Linear Models", arXiv:1712.10024

The replicator should expect the code to run for each table around 10 mins.

The following packages must be installed:

    R 3.6.1 (code was last run with version 3.6.1)
    
    – ```hdm``` (as of 2021-02-16)
    
    – ```xtable``` (as of 2021-02-16)
 
### Replicate code

module load R


Rscript  Run_oracle_std.R Run_oracle_prop.R Run_lasso.R Run_series.R


Rscript print_main_table.R

### Contact

Vira Semenova, semenovavira@gmail.com
