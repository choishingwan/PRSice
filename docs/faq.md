# Frequently Asked Questions
We will continue to update this list to address the more common questions.

1. **I've receive the following error message, what should I do?**

    > *GLM model did not converge! Please send me the DEBUG files*

    This error message means that the logistic regression model cannot
    converge. This is usually caused by small sample size or caused by 
    problem in the input file. 
    
    You should first check the DEBUG file 
    and see if that contains any *NaN* or *Inf*. These will likely be
    caused by un-quality controlled input, which can contain complete
    missingness. If that isn't the case, then you can load the DEBUG and 
    DEBUG.y file into R and see if you can perform the logistic regression
    on the data (DEBUG.y is the y, whereas DEBUG contains the independent 
    variables, including the intercept). 

2. **My base/ target data do not contain the RS ID. Can I still use PRSice?**

    As of version 2.3.x, PRSice now support chromosome ID via the parameter `--chr-id`. `--chr-id` will automatically generate an ID for each of your SNPs based on user input string, some characters are reserved:
    ```
    c = chromosome
    l = coordinates
    a = effective allele
    b = non-effective allele
    ```
    So `--chr-id c:l-ab` will construct SNP ID as `<chromosome>:<coordinate>-<effective allele><non-effective allele>`

    !!! Note 
        
        It is not case sensitive

