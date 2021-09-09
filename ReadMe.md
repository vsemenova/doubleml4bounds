# doubleml4bounds
This R package implements sharp bounds on treatments effects in the presence of selection/nonresponse bias in randomized control trials. It includes basic Lee bounds  [Semenova (2017)](https://arxiv.org/abs/1712.10024). It compares two methods:

-- series-based support function estimator of [Chandrasekhar et al (2012)](https://arxiv.org/abs/1212.5627)

-- lasso-based double machine learning estimator of [Semenova (2017)](https://arxiv.org/abs/1712.10024)

# Packages and Dependencies

The following packages must be installed:

    R 3.6.1 (code was last run with version 3.6.1)
    
    – ```hdm``` (as of 2021-02-16)
    – ```xtable``` (as of 2021-02-16)

# Replication
```
module load R

Rscript Run_series.R Run_lasso.R Run_oracle_std.R Run_oracle_prop.R

Rscript print_main_table.R

```

# Support
Vira Semenova: semenovavira@gmail.com
