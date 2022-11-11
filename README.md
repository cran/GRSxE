
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GRSxE

<!-- badges: start -->
<!-- badges: end -->

GRSxE is a software package for detecting GxE (gene-environment)
interactions using GRS (genetic risk scores). A GRS is constructed on
the data and evaluated for testing an interaction with an environmentalm
exposure while adjusting for potential confounders. The GRS is
constructed using bagging and evaluated performing OOB (out-of-bag)
predictions such that the full data set can be used for both GRS
construction and GxE interaction testing.

## Installation

You can install the released version of GRSxE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GRSxE")
```

## Example

Here is an example of an epidemiological toy data set consisting of some
SNPs, an environmental covariable and a quantitative outcome/phenotype.

``` r
library(GRSxE)
```

### Data generation

``` r
set.seed(101299)
maf <- 0.25
n.snps <- 50
N <- 2000
X <- matrix(sample(0:2, n.snps * N, replace = TRUE,
                   prob = c((1-maf)^2, 1-(1-maf)^2-maf^2, maf^2)), ncol = n.snps)
colnames(X) <- paste("SNP", 1:n.snps, sep="")
E <- rnorm(N, 20, 10)
E[E < 0] <- 0
```

For illustration purposes, an outcome involving a GxE interaction and an
outcome not containing a GxE interaction are constructed and analyzed.

#### Generate outcome with a GxE interaction

``` r
y.GxE <- -0.75 + log(2) * (X[,"SNP1"] != 0) +
  log(4) * E/20 * (X[,"SNP2"] != 0 & X[,"SNP3"] == 0) +
  rnorm(N, 0, 2)
```

#### Generate outcome without a GxE interaction

``` r
y.no.GxE <- -0.75 + log(2) * (X[,"SNP1"] != 0) +
  log(4) * E/20 + log(4) * (X[,"SNP2"] != 0 & X[,"SNP3"] == 0) +
  rnorm(N, 0, 2)
```

### Test for GxE interaction

The GxE test can now be performed by applying the `GRSxE` function.
Since a GLM (generalized linear model) is returned, detailed results can
be retrieved through `summary(...)`.

#### Outcome showing a GxE interaction

First, the outcome involving a GxE interaction is tested.

``` r
summary(GRSxE(X, y.GxE, E))
#> 
#> Call:
#> glm(formula = as.formula(form), family = glm.family, data = dat)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -7.3934  -1.2871  -0.0123   1.3691   6.8683  
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.513301   0.103061  -4.981 6.88e-07 ***
#> G            0.521352   0.280726   1.857   0.0634 .  
#> E            0.028518   0.004626   6.165 8.52e-10 ***
#> G:E          0.055216   0.012654   4.363 1.35e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 3.871335)
#> 
#>     Null deviance: 8574.0  on 1999  degrees of freedom
#> Residual deviance: 7727.2  on 1996  degrees of freedom
#> AIC: 8388.9
#> 
#> Number of Fisher Scoring iterations: 2
```

The corresponding p-value (`G:E`) is very low, indicating there is a GxE
interaction.

#### Outcome not showing a GxE interaction

Next, the outcome not containing a GxE interaction is tested.

``` r
summary(GRSxE(X, y.no.GxE, E))
#> 
#> Call:
#> glm(formula = as.formula(form), family = glm.family, data = dat)
#> 
#> Deviance Residuals: 
#>    Min      1Q  Median      3Q     Max  
#> -7.609  -1.439  -0.022   1.446   6.906  
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -1.9775535  0.3784215  -5.226 1.92e-07 ***
#> G            1.5660013  0.2953620   5.302 1.27e-07 ***
#> E            0.0634424  0.0172752   3.672 0.000247 ***
#> G:E          0.0002192  0.0135055   0.016 0.987054    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 4.453263)
#> 
#>     Null deviance: 10339.6  on 1999  degrees of freedom
#> Residual deviance:  8888.7  on 1996  degrees of freedom
#> AIC: 8669
#> 
#> Number of Fisher Scoring iterations: 2
```

The corresponding p-value (`G:E`) is rather high, leaving no evidence
that there might be a GxE interaction.
