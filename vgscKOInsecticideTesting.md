vgscKOInsecticideTesting
================

``` r
library(lme4); packageVersion('lme4')
```

    ## Loading required package: Matrix

    ## [1] '1.1.38'

``` r
library(ggplot2); packageVersion('ggplot2')
```

    ## [1] '4.0.1'

``` r
library(dplyr); packageVersion('dplyr')
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## [1] '1.1.4'

``` r
library(tidyr); packageVersion('tidyr')
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following objects are masked from 'package:Matrix':
    ## 
    ##     expand, pack, unpack

    ## [1] '1.3.2'

``` r
library(scales); packageVersion('scales')
```

    ## [1] '1.4.0'

``` r
library(ecotox); packageVersion('ecotox')
```

    ## [1] '1.4.4'

``` r
library(arm); packageVersion('arm')
```

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 
    ## arm (Version 1.14-4, built: 2024-4-1)

    ## Working directory is /Users/poppypescod/Library/CloudStorage/OneDrive-LSTM/FUNC GEN GROUP/vgsc gene drive lines/vgscKO/Coding

    ## 
    ## Attaching package: 'arm'

    ## The following object is masked from 'package:scales':
    ## 
    ##     rescale

    ## [1] '1.14.4'

``` r
library(emmeans); packageVersion('emmeans')
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

    ## [1] '2.0.1'

## Deltamethrin bottle assays

``` r
df_Delta <- read.csv("DeltaDeadAlive.csv")
head(df_Delta)
```

    ##   Dose Bottle Genotype Dead Total Alive
    ## 1 0.05      A   vgscKO    1    12    11
    ## 2 0.05      B   vgscKO    0    14    14
    ## 3 0.05      C   vgscKO    0    10    10
    ## 4 0.10      A   vgscKO    0    10    10
    ## 5 0.10      B   vgscKO    0    12    12
    ## 6 0.10      C   vgscKO    0    14    14

GLM with probit link

``` r
Delta <- glmer(cbind(Dead,Alive) ~ log10(Dose) * Genotype + (1 | Dose:Bottle),
              data = df_Delta,
              family = binomial(link = "probit"))

summary(Delta)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( probit )
    ## Formula: cbind(Dead, Alive) ~ log10(Dose) * Genotype + (1 | Dose:Bottle)
    ##    Data: df_Delta
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##      91.2      97.8     -40.6      81.2        23 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.68425 -0.52263 -0.04081  0.69306  2.39536 
    ## 
    ## Random effects:
    ##  Groups      Name        Variance Std.Dev.
    ##  Dose:Bottle (Intercept) 0.1249   0.3534  
    ## Number of obs: 28, groups:  Dose:Bottle, 14
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              0.1602     0.1886   0.849    0.396    
    ## log10(Dose)              2.1405     0.3408   6.280 3.38e-10 ***
    ## GenotypeWT               0.9132     0.2316   3.943 8.04e-05 ***
    ## log10(Dose):GenotypeWT   0.1967     0.4108   0.479    0.632    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) lg10(D) GntyWT
    ## log10(Dose)  0.368               
    ## GenotypeWT  -0.428 -0.172        
    ## lg10(D):GWT -0.122 -0.630   0.538

Check for overdispersion

``` r
rp <- residuals(Delta, type = "pearson")
Delta_over <- sum(rp^2)/df.residual(Delta)
Delta_over
```

    ## [1] 1.01614

Check residuals (should be \<2 for good fit)

``` r
Delta_Model_DevResid <- resid (Delta, type = "deviance")
hist(Delta_Model_DevResid)
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Check variance of residuals across the range of fitted values (most
points should fall within the grey line)

``` r
Delta_Model_Predict <- predict(Delta)
binnedplot(Delta_Model_Predict, Delta_Model_DevResid)
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Check for outliers (Cook’s distance \>1 indicates outliers)

``` r
plot(Delta, which = 4)
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Model fit looks good

### Pairwise comparison of genotype impact

Overall impact of genotype

``` r
emmeans(Delta, pairwise ~ Genotype)
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $emmeans
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO    0.221 0.192 Inf    -0.156     0.598
    ##  WT        1.140 0.234 Inf     0.681     1.599
    ## 
    ## Results are given on the probit (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.919 0.238 Inf  -3.859  0.0001
    ## 
    ## Note: contrasts are still on the probit scale. Consider using
    ##       regrid() if you want contrasts of back-transformed estimates.

Impact of genotype at different insecticide levels

``` r
emmeans(Delta, pairwise ~ Genotype | Dose,
        at = list(Dose = c(0.05, 0.1, 0.5, 1, 5)))
```

    ## $emmeans
    ## Dose = 0.05:
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -2.625 0.413 Inf   -3.4346    -1.815
    ##  WT       -1.967 0.318 Inf   -2.5901    -1.345
    ## 
    ## Dose = 0.10:
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -1.980 0.323 Inf   -2.6139    -1.347
    ##  WT       -1.264 0.238 Inf   -1.7312    -0.797
    ## 
    ## Dose = 0.50:
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -0.484 0.179 Inf   -0.8342    -0.134
    ##  WT        0.370 0.175 Inf    0.0277     0.712
    ## 
    ## Dose = 1.00:
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO    0.160 0.189 Inf   -0.2095     0.530
    ##  WT        1.073 0.228 Inf    0.6272     1.519
    ## 
    ## Dose = 5.00:
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO    1.656 0.354 Inf    0.9624     2.350
    ##  WT        2.707 0.421 Inf    1.8817     3.532
    ## 
    ## Results are given on the probit (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ## Dose = 0.05:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.657 0.454 Inf  -1.448  0.1477
    ## 
    ## Dose = 0.10:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.717 0.346 Inf  -2.068  0.0386
    ## 
    ## Dose = 0.50:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.854 0.195 Inf  -4.375 <0.0001
    ## 
    ## Dose = 1.00:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.913 0.232 Inf  -3.943 <0.0001
    ## 
    ## Dose = 5.00:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -1.051 0.456 Inf  -2.306  0.0211
    ## 
    ## Note: contrasts are still on the probit scale. Consider using
    ##       regrid() if you want contrasts of back-transformed estimates.

**Interpretation** The vgscKO knockout reduces mortality from low dose
deltamethrin exposure significantly overall, and at all doses except
0.05 ug/ml.

### Produce prediction grid for plotting

``` r
df_Delta$BottleID <- paste0(format(df_Delta$Dose, trim = TRUE), df_Delta$Bottle)
head(df_Delta)
```

    ##   Dose Bottle Genotype Dead Total Alive BottleID
    ## 1 0.05      A   vgscKO    1    12    11    0.05A
    ## 2 0.05      B   vgscKO    0    14    14    0.05B
    ## 3 0.05      C   vgscKO    0    10    10    0.05C
    ## 4 0.10      A   vgscKO    0    10    10    0.10A
    ## 5 0.10      B   vgscKO    0    12    12    0.10B
    ## 6 0.10      C   vgscKO    0    14    14    0.10C

``` r
pred_grid <- expand_grid(
    Dose = sort(unique(df_Delta$Dose)),
    Genotype = unique(df_Delta$Genotype)
) |>
    mutate(logdose = log10(Dose))


head(pred_grid)
```

    ## # A tibble: 6 × 3
    ##    Dose Genotype logdose
    ##   <dbl> <chr>      <dbl>
    ## 1  0.05 vgscKO    -1.30 
    ## 2  0.05 WT        -1.30 
    ## 3  0.1  vgscKO    -1    
    ## 4  0.1  WT        -1    
    ## 5  0.5  vgscKO    -0.301
    ## 6  0.5  WT        -0.301

``` r
pred_link <- predict(
  Delta, newdata = pred_grid, type = "link",
  se.fit = TRUE, re.form = NA
)
```

    ## Warning in predict.merMod(Delta, newdata = pred_grid, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
inv_link <- pnorm

pred_df <- pred_grid |>
    mutate(
        fit_link = pred_link$fit,
        se_link = pred_link$se.fit,
        fit = inv_link(fit_link),
        lwr = inv_link(fit_link - 1.96 * se_link),
        upr = inv_link(fit_link + 1.96 * se_link)
        )

# Aggregate observed proportions
obs_df <- df_Delta |>
    group_by(Dose, Genotype) |>
    summarise(Dead = sum(Dead),Total = sum(Total), prop = Dead/Total, .groups = "drop")
```

Plot

``` r
genotype_cols2 <- c(
  "vgscKO" = "#0B3C8C",
  "WT" = "#B11226"
)

genotype_names <- c(
  "vgscKO" = expression(vgscKO^"+"),
  "WT" = expression(vgscKO^"-")
)

ggplot() +
    geom_ribbon(data = pred_df, aes(Dose, ymin = lwr, ymax = upr, fill = Genotype), alpha = 0.18) +
    geom_line(data = pred_df, aes(Dose, y = fit, colour = Genotype), linewidth = 1) +
    geom_point(data = obs_df, aes(Dose, y = prop, colour = Genotype), size = 2.5) +
    scale_color_manual(values = genotype_cols2, labels = genotype_names) +
    scale_fill_manual(values = genotype_cols2, labels = genotype_names) +
    scale_x_log10(breaks = sort(unique(df_Delta$Dose))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = "Deltamethrin dose (µg/ml)", y = "Mortality (%)") +
    theme_classic() +
    theme(legend.position = "bottom")
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## DDT

Load data

``` r
df_DDT <- read.csv("DDTDeadAlive.csv")
head(df_DDT)
```

    ##   Dose Bottle Genotype Dead Total Alive
    ## 1  0.5      A   vgscKO    0     9     9
    ## 2  0.5      B   vgscKO    0     9     9
    ## 3  0.5      C   vgscKO    0     3     3
    ## 4  1.0      A   vgscKO    0     9     9
    ## 5  1.0      B   vgscKO    1    10     9
    ## 6  1.0      C   vgscKO    2     9     7

GLM with probit link

``` r
DDT <- glmer(cbind(Dead,Alive) ~ log10(Dose) * Genotype + (1 | Dose:Bottle),
              data = df_DDT,
              family = binomial(link = "probit"))

summary(DDT)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( probit )
    ## Formula: cbind(Dead, Alive) ~ log10(Dose) * Genotype + (1 | Dose:Bottle)
    ##    Data: df_DDT
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     142.2     150.7     -66.1     132.2        35 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.60070 -0.59338 -0.04063  0.48871  2.82817 
    ## 
    ## Random effects:
    ##  Groups      Name        Variance Std.Dev.
    ##  Dose:Bottle (Intercept) 0.09767  0.3125  
    ## Number of obs: 40, groups:  Dose:Bottle, 20
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)             -1.7930     0.2827  -6.343 2.26e-10 ***
    ## log10(Dose)              1.0133     0.2047   4.951 7.38e-07 ***
    ## GenotypeWT               0.1973     0.3113   0.634    0.526    
    ## log10(Dose):GenotypeWT   0.3475     0.2275   1.527    0.127    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) lg10(D) GntyWT
    ## log10(Dose) -0.877               
    ## GenotypeWT  -0.746  0.658        
    ## lg10(D):GWT  0.651 -0.712  -0.881

Check for overdispersion

``` r
rp <- residuals(DDT, type = "pearson")
DDT_over <- sum(rp^2)/df.residual(DDT)
DDT_over
```

    ## [1] 0.9242944

Check residuals (should be \<2 for good fit)

``` r
DDT_Model_DevResid <- resid (DDT, type = "deviance")
hist(DDT_Model_DevResid)
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Check variance of residuals across the range of fitted values (most
points should fall within the grey line)

``` r
DDT_Model_Predict <- predict(DDT)
binnedplot(DDT_Model_Predict, (DDT_Model_DevResid))
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Check for outliers (Cook’s distance \>1 indicates outliers)

``` r
plot(DDT, which = 4)
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Model fit looks good

### Pairwise comparison of genotype impact

Overall impact of genotype

``` r
emmeans(DDT, pairwise ~ Genotype)
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $emmeans
    ##  Genotype emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -0.319 0.145 Inf    -0.603   -0.0357
    ##  WT        0.384 0.134 Inf     0.121    0.6458
    ## 
    ## Results are given on the probit (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT   -0.703 0.158 Inf  -4.455 <0.0001
    ## 
    ## Note: contrasts are still on the probit scale. Consider using
    ##       regrid() if you want contrasts of back-transformed estimates.

Impact of genotype at different insecticide levels

``` r
emmeans(DDT, pairwise ~ Genotype | Dose,
        at = list(Dose = c(0.5, 1, 5, 10, 25, 50, 100)))
```

    ## $emmeans
    ## Dose =   0.5:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -2.0981 0.338 Inf   -2.7606  -1.43556
    ##  WT       -2.0054 0.256 Inf   -2.5074  -1.50331
    ## 
    ## Dose =   1.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -1.7930 0.283 Inf   -2.3471  -1.23897
    ##  WT       -1.5957 0.213 Inf   -2.0136  -1.17780
    ## 
    ## Dose =   5.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -1.0848 0.172 Inf   -1.4211  -0.74845
    ##  WT       -0.6445 0.133 Inf   -0.9059  -0.38315
    ## 
    ## Dose =  10.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -0.7797 0.143 Inf   -1.0591  -0.50036
    ##  WT       -0.2349 0.119 Inf   -0.4676  -0.00221
    ## 
    ## Dose =  25.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -0.3765 0.141 Inf   -0.6531  -0.09999
    ##  WT        0.3066 0.130 Inf    0.0525   0.56079
    ## 
    ## Dose =  50.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO   -0.0715 0.169 Inf   -0.4018   0.25885
    ##  WT        0.7163 0.157 Inf    0.4087   1.02385
    ## 
    ## Dose = 100.0:
    ##  Genotype  emmean    SE  df asymp.LCL asymp.UCL
    ##  vgscKO    0.2335 0.211 Inf   -0.1799   0.64699
    ##  WT        1.1259 0.193 Inf    0.7469   1.50492
    ## 
    ## Results are given on the probit (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ## Dose =   0.5:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.0927 0.373 Inf  -0.249  0.8037
    ## 
    ## Dose =   1.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.1973 0.311 Inf  -0.634  0.5262
    ## 
    ## Dose =   5.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.4402 0.187 Inf  -2.355  0.0185
    ## 
    ## Dose =  10.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.5448 0.154 Inf  -3.527  0.0004
    ## 
    ## Dose =  25.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.6831 0.154 Inf  -4.448 <0.0001
    ## 
    ## Dose =  50.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.7878 0.185 Inf  -4.255 <0.0001
    ## 
    ## Dose = 100.0:
    ##  contrast    estimate    SE  df z.ratio p.value
    ##  vgscKO - WT  -0.8924 0.233 Inf  -3.828  0.0001
    ## 
    ## Note: contrasts are still on the probit scale. Consider using
    ##       regrid() if you want contrasts of back-transformed estimates.

**Interpretation** The vgscKO knockout reduces mortality from low dose
DDT exposure significantly overall, and at all doses except 0.5 ug/ml
and 1 ug/ml.

### Produce prediction grid for plotting

``` r
df_DDT$BottleID <- paste0(format(df_DDT$Dose, trim = TRUE), df_DDT$Bottle)
head(df_DDT)
```

    ##   Dose Bottle Genotype Dead Total Alive BottleID
    ## 1  0.5      A   vgscKO    0     9     9     0.5A
    ## 2  0.5      B   vgscKO    0     9     9     0.5B
    ## 3  0.5      C   vgscKO    0     3     3     0.5C
    ## 4  1.0      A   vgscKO    0     9     9     1.0A
    ## 5  1.0      B   vgscKO    1    10     9     1.0B
    ## 6  1.0      C   vgscKO    2     9     7     1.0C

``` r
pred_grid <- expand_grid(
    Dose = sort(unique(df_DDT$Dose)),
    Genotype = unique(df_DDT$Genotype)
) |>
    mutate(logdose = log10(Dose))


head(pred_grid)
```

    ## # A tibble: 6 × 3
    ##    Dose Genotype logdose
    ##   <dbl> <chr>      <dbl>
    ## 1   0.5 vgscKO    -0.301
    ## 2   0.5 WT        -0.301
    ## 3   1   vgscKO     0    
    ## 4   1   WT         0    
    ## 5   5   vgscKO     0.699
    ## 6   5   WT         0.699

``` r
pred_link <- predict(
  DDT, newdata = pred_grid, type = "link",
  se.fit = TRUE, re.form = NA
)
```

    ## Warning in predict.merMod(DDT, newdata = pred_grid, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
inv_link <- pnorm

pred_df <- pred_grid |>
    mutate(
        fit_link = pred_link$fit,
        se_link = pred_link$se.fit,
        fit = inv_link(fit_link),
        lwr = inv_link(fit_link - 1.96 * se_link),
        upr = inv_link(fit_link + 1.96 * se_link)
        )

# Aggregate observed proportions
obs_df <- df_DDT |>
    group_by(Dose, Genotype) |>
    summarise(Dead = sum(Dead),Total = sum(Total), prop = Dead/Total, .groups = "drop")
```

Plot

``` r
genotype_cols2 <- c(
  "vgscKO" = "#0B3C8C",
  "WT" = "#B11226"
)

genotype_names <- c(
  "vgscKO" = expression(vgscKO^"+"),
  "WT" = expression(vgscKO^"-")
)

ggplot() +
    geom_ribbon(data = pred_df, aes(Dose, ymin = lwr, ymax = upr, fill = Genotype), alpha = 0.18) +
    geom_line(data = pred_df, aes(Dose, y = fit, colour = Genotype), linewidth = 1) +
    geom_point(data = obs_df, aes(Dose, y = prop, colour = Genotype), size = 2.5) +
    scale_color_manual(values = genotype_cols2, labels = genotype_names) +
    scale_fill_manual(values = genotype_cols2, labels = genotype_names) +
    scale_x_log10(breaks = sort(unique(df_DDT$Dose))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = "DDT dose (µg/ml)", y = "Mortality (%)") +
    theme_classic() +
    theme(legend.position = "bottom")
```

![](vgscKOInsecticideTesting_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## LC tables for both insecticides

Load combined data

``` r
df <- read.csv("InsecticideDeadAlive.csv")
head(df)
```

    ##    Insecticide Dose Bottle Genotype Dead Total Alive
    ## 1 Deltamethrin 0.05      A   vgscKO    1    12    11
    ## 2 Deltamethrin 0.05      B   vgscKO    0    14    14
    ## 3 Deltamethrin 0.05      C   vgscKO    0    10    10
    ## 4 Deltamethrin 0.10      A   vgscKO    0    10    10
    ## 5 Deltamethrin 0.10      B   vgscKO    0    12    12
    ## 6 Deltamethrin 0.10      C   vgscKO    0    14    14

``` r
df_clean <- subset(df,
                   Dose > 0 &
                   !is.na(Dose) &
                   !is.na(Dead) &
                   !is.na(Total) &
                   Total > 0
)
nrow(df)
```

    ## [1] 72

``` r
nrow(df_clean)
```

    ## [1] 68

``` r
KO_delta <- LC_probit((Dead/Total) ~ log10(Dose),
          p = c(50,90,95,99),
          weights = Total,
          data = df_clean,
          subset = df_clean$Genotype == "vgscKO" & df_clean$Insecticide == "Deltamethrin"
)

WT_delta <- LC_probit((Dead/Total) ~ log10(Dose),
          p = c(50,90,95,99),
          weights = Total,
          data = df_clean,
          subset = df_clean$Genotype == "WT" & df_clean$Insecticide == "Deltamethrin"
)

KO_DDT<- LC_probit((Dead/Total) ~ log10(Dose),
          p = c(50,90,95,99),
          weights = Total,
          data = df_clean,
          subset = df_clean$Genotype == "vgscKO" & df_clean$Insecticide == "DDT"
)

WT_DDT<- LC_probit((Dead/Total) ~ log10(Dose),
          p = c(50,90,95,99),
          weights = Total,
          data = df_clean,
          subset = df_clean$Genotype == "WT" & df_clean$Insecticide == "DDT"
)
```

``` r
KO_delta2 <- KO_delta %>%
  mutate(Genotype = "vgscKO",
         Insecticide = "Deltamethrin")

WT_delta2 <- WT_delta %>%
  mutate(Genotype = "WT",
         Insecticide = "Deltamethrin")

KO_DDT2 <- KO_DDT %>%
  mutate(Genotype = "vgscKO",
         Insecticide = "DDT")

WT_DDT2 <- WT_DDT %>%
  mutate(Genotype = "WT",
         Insecticide = "DDT")

LC_all <- bind_rows(KO_delta2, WT_delta2, KO_DDT2, WT_DDT2)

LC_trimmed <- LC_all %>%
  dplyr::select(Genotype, Insecticide, p, dose, LCL, UCL, se)

LC_trimmed
```

    ## # A tibble: 16 × 7
    ##    Genotype Insecticide      p      dose      LCL        UCL    se
    ##    <chr>    <chr>        <dbl>     <dbl>    <dbl>      <dbl> <dbl>
    ##  1 vgscKO   Deltamethrin    50     0.770    0.455      1.48   1.16
    ##  2 vgscKO   Deltamethrin    90     3.43     1.70      19.0    1.32
    ##  3 vgscKO   Deltamethrin    95     5.24     2.33      41.7    1.39
    ##  4 vgscKO   Deltamethrin    99    11.6      4.11     185.     1.53
    ##  5 WT       Deltamethrin    50     0.351    0.272      0.453  1.14
    ##  6 WT       Deltamethrin    90     1.34     0.953      2.17   1.22
    ##  7 WT       Deltamethrin    95     1.95     1.32       3.48   1.27
    ##  8 WT       Deltamethrin    99     3.97     2.40       8.56   1.36
    ##  9 vgscKO   DDT             50    58.3     34.1      142.     1.39
    ## 10 vgscKO   DDT             90  1103.     342.     12789.     2.23
    ## 11 vgscKO   DDT             95  2538.     636.     47308.     2.60
    ## 12 vgscKO   DDT             99 12121.    2025.    555183.     3.46
    ## 13 WT       DDT             50    15.5      8.59      28.3    1.17
    ## 14 WT       DDT             90   140.      64.3      638.     1.32
    ## 15 WT       DDT             95   261.     104.      1679.     1.40
    ## 16 WT       DDT             99   839.     253.     10586.     1.57
