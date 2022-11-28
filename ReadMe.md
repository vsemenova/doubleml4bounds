# doubleml4bounds
This is replication code for  [Semenova (2017)](https://arxiv.org/abs/1712.10024). It estimates the upper bound on a treatment effect parameter when the outcome variable is recorded in brackets. The code compares two methods:

-- series-based support function estimator of [Chandrasekhar et al (2012)](https://arxiv.org/abs/1212.5627)

-- lasso-based double machine learning estimator of [Semenova (2017)](https://arxiv.org/abs/1712.10024)

# Packages and Dependencies

The following packages must be installed:

    R 3.6.1 (code was last run with version 3.6.1)
    
    – ```hdm``` (as of 2021-02-16)
    – ```xtable``` (as of 2021-02-16)

# Replication of simulations
```
module load R

cd simulations 

Rscript Run_series.R Run_lasso.R Run_oracle_std.R Run_oracle_prop.R

Rscript print_main_table.R

```

# Replication of empirical exercise
```
module load R

cd empirical_applications 

Rscript Run_lasso.R Run_forest.R

Rscript print_main_table.R

```

# Support
Vira Semenova: semenovavira@gmail.com
