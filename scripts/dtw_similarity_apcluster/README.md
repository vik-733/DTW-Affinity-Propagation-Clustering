 
 # Description

 The **`dtw_similarity`** function is a wraper function that enables the use of dynamic time wraping (DTW) distance,
 as implemented in the **`dtwclust`** package, in the context of affinity propagation clustering performed by the **`apcluster`** package.
 The end result is clustering of genomic features based on spatial point sequences along their lengths—such as 100bp bins—in which the values might represent sequencing signal(s).

## Purpose and Scope 

The function is intended for grouping genomic features into clusters based on their spatial profiles, which could include single-variable cases or multivariable scenarios—such as clustering based on profiles across different samples.

## Usage

Once you have sourced the code, clustering can be performed by invoking the **`apcluster`** (or **`apclusterL`**) function and using **`dtw_similarity`** as the 's' argument. Refer to examples in examples.R for more details.

```
library(apcluster)
source(dtw_similarity.R)

cluster_results <- apcluster(s=dtw_similarity,x=data)
```

## Data Format

**1. Univariable Cases:**

- **Option 1** 
    - **List:** One item per genomic feature, vector of values containing the sequence of data points measured along the genomic feature. 

- **Option 2** (Only for same number of bins for all genomic features)

    - **Matrix / Data Frame:** Rows should correspond to individual genomic features, with columns containing the sequence of data points measured along the genomic feature. 

**2. Multivariable cases:**

- **List of Matrices:** When dealing with data from multiple variables, this structure is needed:
  - Each matrix within the list should cointain measurements for an individual genomic feature. 
  - Rows in the matrix are the ordered sequence of spatial intervals or positional bins along the genomic feature. 
  - Each column should represent one of the variables. 

See also examples in in examples.R.

