# jafar

**Bayesian Joint Additive Factor Regression for Multi-View Learning**

`jafar` is an `R` package implementing supervised Bayesian factor models for multi-view data integration.
The associated metodologies, introduced in [Anceschi et al. (2025)](https://arxiv.org/abs/2406.00778), include 

- **Joint Factor Regression ($\texttt{JFR}$):** a baseline model that captures the combined variation across multiple data views using a single set of latent factors.
- **Joint Additive Factor Regression ($\texttt{JAFAR}$):** a more refined model that explicitly decomposes variation into shared and view-specific components.

Both models leverage extensions of the $\texttt{CUSP}$ prior [(Legramanti et al., 2020)](http://academic.oup.com/biomet/article/107/3/745/5847840), learning adaptively the number of factors in a fully Bayesian way.

The package supports posterior sampling for continuous predictors and both continuous and binary response.  
It includes data preprocessing, results visualization, and post-processing functions (e.g. for resolving rotational ambiguity).  

## Installation

```r
# From GitHub - requires devtools
devtools::install_github("niccoloanceschi/jafar")
```
## Documentation and Tutorial

The main functions are detailed in the package [manual](https://github.com/niccoloanceschi/jafar/blob/jafar-manual.pdf).

A full [tutorial](https://github.com/niccoloanceschi/jafar/blob/main/tutorial/tutorial.html) with examples and detailed workflow is provided in the repository.

## Example usage 

```{r}
library(jafar)

# Load data
data <- readRDS("path/to/data.rds")

# Preprocess
X_m <- preprocess_X(data$X_m)
y <- preprocess_y(data$yTrain)

# Fit JAFAR
fit <- gibbs_jafar(X$X_m, y = y$yTrain,
                   K0 = 25, K0_m = rep(20, length(X_m)),
                   tMCMC = 1000, tBurnIn = 500, tThin = 1)

# Predictions
y_pred <- predict_y(X_m, fit)
```

## Reference

Niccolo Anceschi, Federico Ferrari, David B. Dunson, & Himel Mallick. (2025). Bayesian Joint Additive Factor Models for Multiview Learning. [arXiv:2406.00778](https://arxiv.org/abs/2406.00778), include 


