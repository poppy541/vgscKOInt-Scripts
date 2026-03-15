vgscKO fecundity analysis
================

``` r
library(emmeans); packageVersion("emmeans")
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

    ## [1] '2.0.1'

``` r
library(glmmTMB); packageVersion("glmmTMB")
```

    ## [1] '1.1.14'

``` r
library(MASS); packageVersion("MASS")
```

    ## [1] '7.3.65'

``` r
library(arm); packageVersion("arm")
```

    ## Loading required package: Matrix

    ## Loading required package: lme4

    ## 
    ## arm (Version 1.14-4, built: 2024-4-1)

    ## Working directory is /Users/poppypescod/Library/CloudStorage/OneDrive-LSTM/FUNC GEN GROUP/vgsc gene drive lines/vgscKO/Coding

    ## [1] '1.14.4'

``` r
library(DHARMa); packageVersion("DHARMa")
```

    ## This is DHARMa 0.4.7. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')

    ## [1] '0.4.7'

Loading data

``` r
df <- read.csv("HatchData.csv")
head(df)
```

    ##   vgscKOParent Eggs Larvae Hatch Genotype
    ## 1         male  140    116  82.9   vgscKO
    ## 2         male   99     93  93.9       WT
    ## 3         male   96     20  20.8   vgscKO
    ## 4         male  137     90  65.7   vgscKO
    ## 5         male   82     77  93.9       WT
    ## 6         male   72     34  47.2       WT

Removing batches with fewer than 10 larvae

``` r
df_clean <- na.omit(df)
df_clean <- subset(df_clean, Larvae > 10)

df_clean$Genotype <- factor(df_clean$Genotype, levels = c("WT", "vgscKO"))
df_clean$vgscKOParent <- factor(df_clean$vgscKOParent, levels = c("female", "male"))
```

### Egg laying

Negative binomial GLM to determine impact of Genotype and vgscKOParent
on number of eggs laid

``` r
m_eggs <- glm.nb(Eggs ~ Genotype * vgscKOParent, data = df_clean)
summary(m_eggs)
```

    ## 
    ## Call:
    ## glm.nb(formula = Eggs ~ Genotype * vgscKOParent, data = df_clean, 
    ##     init.theta = 10.73838928, link = log)
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                      4.89622    0.06614  74.034  < 2e-16 ***
    ## GenotypevgscKO                  -0.18642    0.08589  -2.170 0.029978 *  
    ## vgscKOParentmale                -0.36549    0.09429  -3.876 0.000106 ***
    ## GenotypevgscKO:vgscKOParentmale  0.32852    0.12584   2.611 0.009038 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(10.7384) family taken to be 1)
    ## 
    ##     Null deviance: 124.41  on 105  degrees of freedom
    ## Residual deviance: 108.92  on 102  degrees of freedom
    ## AIC: 1060.4
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  10.74 
    ##           Std. Err.:  1.62 
    ## 
    ##  2 x log-likelihood:  -1050.37

``` r
m_eggs_no_int <- glm.nb(Eggs ~ Genotype + vgscKOParent, data = df_clean)
anova(m_eggs_no_int, m_eggs, test = "Chisq")
```

    ## Likelihood ratio tests of Negative Binomial Models
    ## 
    ## Response: Eggs
    ##                     Model    theta Resid. df    2 x log-lik.   Test    df
    ## 1 Genotype + vgscKOParent 10.02762       103       -1056.984             
    ## 2 Genotype * vgscKOParent 10.73839       102       -1050.370 1 vs 2     1
    ##   LR stat.    Pr(Chi)
    ## 1                    
    ## 2 6.613771 0.01011932

Model with interaction term kept

Check for overdispersion

``` r
dispersion_eggs <- sum(residuals(m_eggs, type = "pearson")^2) / df.residual(m_eggs)
dispersion_eggs
```

    ## [1] 0.8697583

Check residuals (should be \<2 for good fit)

``` r
m_eggs_Model_DevResid <- resid (m_eggs, type = "deviance")
hist(m_eggs_Model_DevResid)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Check variance of residuals across the range of fitted values (most
points should fall within the grey line)

``` r
m_eggs_Model_Predict <- predict(m_eggs)
binnedplot(m_eggs_Model_Predict, (m_eggs_Model_DevResid))
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Check for outliers (Cook’s distance \>1 indicates outliers)

``` r
plot(m_eggs, which = 4)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Impact of genotype on egg production

``` r
pairs(emmeans(m_eggs, ~ Genotype, weights = "proportional"))
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ##  contrast    estimate     SE  df z.ratio p.value
    ##  WT - vgscKO   0.0346 0.0628 Inf   0.550  0.5820
    ## 
    ## Results are averaged over the levels of: vgscKOParent 
    ## Results are given on the log (not the response) scale.

**Interpretation**

vgscKO females produce slightly larger egg batches than vgscKO males -
but there is no overall impact of the vgscKO knock out on egg
production.

### Larvae numbers

Negative binomial GLM for larvae numbers, with an offset factor for egg
number.

``` r
m_larvae_eff <- glm.nb(
  Larvae ~ Genotype * vgscKOParent + offset(log(Eggs)),
  data = df_clean
)
summary(m_larvae_eff)
```

    ## 
    ## Call:
    ## glm.nb(formula = Larvae ~ Genotype * vgscKOParent + offset(log(Eggs)), 
    ##     data = df_clean, init.theta = 11.99897303, link = log)
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     -0.32356    0.06407  -5.050 4.43e-07 ***
    ## GenotypevgscKO                  -0.02445    0.08360  -0.292    0.770    
    ## vgscKOParentmale                -0.05265    0.09209  -0.572    0.567    
    ## GenotypevgscKO:vgscKOParentmale  0.07732    0.12288   0.629    0.529    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(11.999) family taken to be 1)
    ## 
    ##     Null deviance: 109.22  on 105  degrees of freedom
    ## Residual deviance: 108.76  on 102  degrees of freedom
    ## AIC: 972.03
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  12.00 
    ##           Std. Err.:  1.94 
    ## 
    ##  2 x log-likelihood:  -962.026

``` r
dispersion_larvae <- sum(residuals(m_larvae_eff, type = "pearson")^2) / df.residual(m_larvae_eff)
dispersion_larvae
```

    ## [1] 0.8142196

Check residuals (should be \<2 for good fit)

``` r
larvae_Model_DevResid <- resid (m_larvae_eff, type = "deviance")
hist(larvae_Model_DevResid)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Slightly skewed but still relatively normal distribution around 0

Check variance of residuals across the range of fitted values (most
points should fall within the grey line)

``` r
larvae_Model_Predict <- predict(m_larvae_eff)
binnedplot(larvae_Model_Predict, (larvae_Model_DevResid))
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Check for outliers (Cook’s distance \>1 indicates outliers)

``` r
plot(m_larvae_eff, which = 4)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Double-check with DHARMa

``` r
sim <- simulateResiduals(m_larvae_eff)
plot(sim)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
testDispersion(sim)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.83512, p-value = 0.424
    ## alternative hypothesis: two.sided

Model fits well

Estimated marginal means

``` r
pairs(emmeans(m_larvae_eff, ~ Genotype, type = "response"))
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ##  contrast    ratio     SE  df null z.ratio p.value
    ##  WT / vgscKO 0.986 0.0606 Inf    1  -0.231  0.8171
    ## 
    ## Results are averaged over the levels of: vgscKOParent 
    ## Tests are performed on the log scale

``` r
pairs(emmeans(m_larvae_eff, ~ vgscKOParent, type = "response"))
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ##  contrast      ratio     SE  df null z.ratio p.value
    ##  female / male  1.01 0.0623 Inf    1   0.228  0.8199
    ## 
    ## Results are averaged over the levels of: Genotype 
    ## Tests are performed on the log scale

**Interpretation**

No significant impact of genotype or vgscKO parent on larval output.

### Hatch rate analysis

Quasibinomial GLM

``` r
df_hatch <- subset(df_clean, Eggs > 0 & Larvae <= Eggs)

m_hatch_q <- glm(
  cbind(Larvae, Eggs - Larvae) ~ Genotype * vgscKOParent,
  family = quasibinomial,
  data = df_hatch
)
summary(m_hatch_q)
```

    ## 
    ## Call:
    ## glm(formula = cbind(Larvae, Eggs - Larvae) ~ Genotype * vgscKOParent, 
    ##     family = quasibinomial, data = df_hatch)
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                      0.93410    0.18747   4.983 2.57e-06 ***
    ## GenotypevgscKO                  -0.06463    0.25104  -0.257    0.797    
    ## vgscKOParentmale                -0.08122    0.29000  -0.280    0.780    
    ## GenotypevgscKO:vgscKOParentmale  0.21577    0.38998   0.553    0.581    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 21.90005)
    ## 
    ##     Null deviance: 2383.8  on 105  degrees of freedom
    ## Residual deviance: 2375.9  on 102  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
dispersion_hatch <- sum(residuals(m_hatch_q, type = "pearson")^2) / df.residual(m_hatch_q)
dispersion_hatch
```

    ## [1] 21.89997

Dispersion far too high - try a beta-binomial model

``` r
library(glmmTMB)
m_hatch_bb <- glmmTMB(
  cbind(Larvae, Eggs - Larvae) ~ Genotype * vgscKOParent,
  family = betabinomial(link = "logit"),
  data = df_hatch
)
summary(m_hatch_bb)
```

    ##  Family: betabinomial  ( logit )
    ## Formula:          cbind(Larvae, Eggs - Larvae) ~ Genotype * vgscKOParent
    ## Data: df_hatch
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     926.9     940.2    -458.4     916.9       101 
    ## 
    ## 
    ## Dispersion parameter for betabinomial family (): 4.28 
    ## 
    ## Conditional model:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                       1.0789     0.2028   5.320 1.04e-07 ***
    ## GenotypevgscKO                   -0.2041     0.2567  -0.795    0.427    
    ## vgscKOParentmale                 -0.3062     0.2772  -1.104    0.269    
    ## GenotypevgscKO:vgscKOParentmale   0.3301     0.3676   0.898    0.369    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

**Check for model fit using DHARMa**

Overdispersion

``` r
sim <- simulateResiduals(m_hatch_bb)
plot(sim)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
testDispersion(sim)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.95548, p-value = 0.688
    ## alternative hypothesis: two.sided

Check for outliers

``` r
testOutliers(sim)
```

![](vgscKOanalysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

    ## 
    ##  DHARMa bootstrapped outlier test
    ## 
    ## data:  sim
    ## outliers at both margin(s) = 0, observations = 106, p-value = 0.8
    ## alternative hypothesis: two.sided
    ##  percent confidence interval:
    ##  0.00000000 0.02830189
    ## sample estimates:
    ## outlier frequency (expected: 0.00792452830188679 ) 
    ##                                                  0

Beta-binomial model fits well

Overall (marginal) genotype effect:

``` r
pairs(emmeans(m_hatch_bb, ~ Genotype, type = "response"))
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ##  contrast    odds.ratio    SE  df null z.ratio p.value
    ##  WT / vgscKO       1.04 0.191 Inf    1   0.212  0.8318
    ## 
    ## Results are averaged over the levels of: vgscKOParent 
    ## Tests are performed on the log odds ratio scale

``` r
pairs(emmeans(m_hatch_bb, ~ vgscKOParent, type = "response"))
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ##  contrast      odds.ratio    SE  df null z.ratio p.value
    ##  female / male       1.15 0.212 Inf    1   0.768  0.4424
    ## 
    ## Results are averaged over the levels of: Genotype 
    ## Tests are performed on the log odds ratio scale

**Interpretation**

There was no significant effect of Genotype or vgscKO Parent sex on the
hatch rate of eggs.

# Likelihood analysis for homozygotes in group lays

``` r
# Data
pos <- c(338, 513, 448)
tot <- c(500, 759, 661)

# Hypotheses
pA <- 0.75 # Proportion of clutch positive if homs exist
pB <- 0.67 # Proportion of clutch positive if homs don't exist
```

``` r
# Log-likelihood under binomial for each hypothesis
llA <- sum(dbinom(pos, size = tot, prob = pA, log = TRUE))
llB <- sum(dbinom(pos, size = tot, prob = pB, log = TRUE))

# Likelihood ratio in favour of B over A
LR_B_over_A <- exp(llB - llA)

# Report in log10 (Jeffreys-style scale)
log10_LR_B_over_A <- (llB - llA) / log(10)
list(
  logLik_A = llA,
  logLik_B = llB,
  LR_B_over_A = LR_B_over_A,
  log10_LR_B_over_A = log10_LR_B_over_A
)
```

    ## $logLik_A
    ## [1] -36.22815
    ## 
    ## $logLik_B
    ## [1] -10.34019
    ## 
    ## $LR_B_over_A
    ## [1] 174984202038
    ## 
    ## $log10_LR_B_over_A
    ## [1] 11.243

**Interpretation**

The data are 1.75 x 10^11 more likely to occur with a
67% ratio (e.g. homozygotes don’t exist) than a 75% ratio (they do
exist).
