# Empirical P-value calculation
All approaches to PRS calculation involve parameter optimisation in generating the final prediction model, and are thus vulnerable to overfitting.
There are a few methods to account for the overfitting:

1. Evaluate performance in an independent validation sample
2. Cross validation
3. Calculate an empirical P-value

In, PRSice-2, we have implemented a permutation analysis to calculate an empirical P-value. 

# Permutation Procedure
To calculate the empirical P-value, PRSice-2 perform the following

1. Perform standard PRSice analysis 
    - Obtain the p-value of association of the best p-value threshold ($P_o$)
2. Randomly shuffle the phenotype and repeat the PRSice analysis 
    - Obtain the p-value of association of the best p-value threshold under the null ($P_{null}$)
3. Repeat step-2 $N=10,000$ times (for `--perm 10000`)
4. The empirical p-value can then be calculated as

$$
\text{Empirical-}P = \frac{\sum_{n=1}^NI(P_{null}\lt P_o)+1}{N+1}
$$

where $I(.)$ is the indicator function. 

!!! Warning
    
    While the empirical p-value for association will be controlled for Type 1 error, 
    the observed phenotypic variance explained, R^2^, remains unadjusted and is affected by overfitting. 
    Therefore, it is imperative to perform out-of-samp,le prediction, or cross-validation to evaluate the predictive accuracy of PRS. 

# Computation Algorithm
In reality, due to the time consuming nature of permutation analysis, PRSice-2 exploit certain property of random number generation to speed up the permutation analysis. 

To generate random numbers, a random seed is required. When the same seed is provided, the same sequence of random number will always be generated. 
PRSice-2 exploit this property, such that the permutation analysis is performed as follow

1. Generate the random seed or use the user provided random seed to $S$
2. For each p-value threshold
    1. Calculate the observed p-value 
    2. Seed the random number generator with $S$
    3. For Quantitative trait, (and binary trait, unless `--logit-perm` is set), decompose the matrix of the independent variables ($Intercept+PRS+Covariates$)
    4. Generate N copies of random phenotypes via random shuffling. 
    5. Calculate the p-value association for each null phenotype
    6. For each permutation, check if the current null p-value is the most significant. Replace the previous "best" p-value if the current null p-value is more significant
3. Calculate the empirical p-value once all p-value thresholds have been processed

As we re-seed the random number generator for each p-value threshold, we ensure the random phenotypes generated in each p-value thresholds
to be identical, therefore allowing us to reuse the calculated PRS, and also reuse the decomosed matrix, which leads to significant speed 
up of the permutation process.