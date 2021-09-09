The code in this replication package replicates the simulation results contained in Tables 1,2 for the paper

Semenova, Vira "Debiased Machine Learning of Set-Identified Linear Models", arXiv:1712.10024

The replicator should expect the code to run for each table around 10 mins.


The following packages must be installed:
    R 3.6.1 (code was last run with version 3.6.1)
    – hdm (as of 2021-02-16)
    – xtable (as of 2021-02-16)

Tables from Main Text and Appendix are referred to as MainTables. The intermediate tables are saved as Tablek.csv, for  k in {oracle_std, oracle_prop, lasso, series}.  


1. Open  Runk.R, where  k in {oracle_std, oracle_prop, lasso, series}.  

2. Adjust the foldername in line 6

3. Run Runk.R for each k. The resulting Tablek.csv is saved in ../Tables

4. Run print_main_table.R to get MainTable.txt (Table 1, main text) and AppendixTable.txt (Table A.1, appendix).


Step 3 details:

The output of Runk.R is the risk and the rejection frequency of the output of the Algorithm 1. The Tables differ in terms of the input values of the algorithm.


Results:
    Table_oracle_std.csv:  true (treatment fitted values. zero outcome fitted values (i.e., partialling-out=FALSE)  
    
    Table_oracle_prop.csv: true treatment fitted values. LASSO outcome fitted values (i.e., partialling-out=TRUE)  
             
    Table_lasso.csv: LASSO treatment fitted values. LASSO outcome fitted values (i.e., partialling-out=TRUE). performance close to oracle (Table 2). 
             
    Table_series.csv series treatment fitted values. zero outcome fitted values (i.e., partialling-out=FALSE). Coincides with Chandrasekhar et al (2012).        



Important functions:

Functions.R

    bracket -- assign an observed outcome into brackets of width Delta

    simulate_data_2d -- generate data (the 2-dimensional treatment vector, the outcome, and its brackets) as described in Section 5 of the paper

    estimate_upper_bound -- encodes the Algorithm 1 (1-dimensional case) of the paper given the first-stage fitted values

    estimate_support_function -- encodes the Algorithm 1 (2-dimensional case) of the paper given the first-stage fitted values

    
FirstStage.R
    firststage.R -- predict the treatment vector and the lower bound bracket given method_treat, method_outcome using cross-fitting with K=2 folds.

RemoveControls.R
   contains a list of helper functions that encode cross-fitting
   

    
