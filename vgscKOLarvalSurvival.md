vgscKO larval survival
================

``` r
library(dplyr); packageVersion("dplyr")
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
library(lme4); packageVersion("lme4")
```

    ## Loading required package: Matrix

    ## [1] '1.1.38'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '4.0.1'

``` r
library(arm); packageVersion("arm")
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

    ## [1] '1.14.4'

``` r
library(emmeans); packageVersion("emmeans")
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

    ## [1] '2.0.1'

``` r
library(DHARMa); packageVersion("DHARMa")
```

    ## This is DHARMa 0.4.7. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')

    ## [1] '0.4.7'

Load data

``` r
df <- read.csv("SurvivalRate.csv")
head(df)
```

    ##   Genotype Tray NoLarveStart Day Larvae Perc.survival
    ## 1     CFP+    A          25L   1     25           100
    ## 2     CFP+    B          25L   1     25           100
    ## 3     CFP+    C          25L   1     25           100
    ## 4     CFP+    D          25L   1     25           100
    ## 5     CFP+    E          50L   1     50           100
    ## 6     CFP+    F          50L   1     50           100

## Binomial GLMM with complementary log-log link

First reformat data to compute at risk numbers (necessary for the GLMM)

``` r
cp <- df %>%
  # rename columns
  rename(
    genotype  = Genotype,
    tray_raw  = Tray,
    tray_size = NoLarveStart,
    day       = Day,
    alive     = Larvae
  ) %>%
  
  # Create unique tray IDs (since A, B, C appear in each genotype)
  mutate(
    genotype  = factor(genotype),
    tray_size = factor(tray_size, levels = c("25L", "50L")),
    day       = as.integer(day),
    tray_id   = interaction(genotype, tray_raw, drop = TRUE, sep = "_")
  ) %>%
  
  arrange(genotype, tray_id, day) %>%
  group_by(genotype, tray_id, tray_size) %>%
  
  # Compute alive_next and deaths
  mutate(
    alive_next = lead(alive),

    # For the last observed day:
    # larvae are censored (so alive_next = alive → deaths = 0)
    alive_next = ifelse(is.na(alive_next), alive, alive_next),

    # Number dying in the interval (day to day+1)
    deaths = pmax(alive - alive_next, 0L),

    # Number at risk at the start of the day
    at_risk = alive
  ) %>%
  ungroup()

cp <- cp %>%
  mutate(genotype = relevel(genotype, ref = "G3"))  # sets G3 as control
```

Binomial GLMM

``` r
m1 <- glmer(
  cbind(deaths, at_risk - deaths) ~ genotype + tray_size + factor(day) + (1 | tray_id),
  data = cp,
  family = binomial(link = "cloglog"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
summary(m1)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( cloglog )
    ## Formula: 
    ## cbind(deaths, at_risk - deaths) ~ genotype + tray_size + factor(day) +  
    ##     (1 | tray_id)
    ##    Data: cp
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     285.0     329.7    -128.5     257.0       166 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5322 -0.5925 -0.3521  0.0000  3.7601 
    ## 
    ## Random effects:
    ##  Groups  Name        Variance Std.Dev.
    ##  tray_id (Intercept) 0        0       
    ## Number of obs: 180, groups:  tray_id, 18
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -3.2790     0.4083  -8.031 9.64e-16 ***
    ## genotypeCFP-     0.1487     0.3690   0.403 0.686985    
    ## genotypeCFP+    -0.2264     0.3838  -0.590 0.555299    
    ## tray_size50L     0.1974     0.2484   0.795 0.426893    
    ## factor(day)2    -0.7368     0.3358  -2.194 0.028208 *  
    ## factor(day)3    -1.3398     0.4226  -3.170 0.001522 ** 
    ## factor(day)4    -1.6657     0.4855  -3.431 0.000601 ***
    ## factor(day)5   -20.9237  6812.9227  -0.003 0.997550    
    ## factor(day)6    -2.5758     0.7320  -3.519 0.000433 ***
    ## factor(day)7    -1.3149     0.4226  -3.112 0.001859 ** 
    ## factor(day)8    -0.9422     0.3683  -2.558 0.010531 *  
    ## factor(day)9    -1.4381     0.4499  -3.196 0.001392 ** 
    ## factor(day)10  -20.8856  6839.8813  -0.003 0.997564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 13 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

    ## optimizer (bobyqa) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

Model checks - whether genotype impacts the model, whether there are
interactions between tray and genotype, and overdispersion checks

``` r
# Comparing models with and without genotype (to see if genotype makes a difference)
m0 <- update(m1, . ~ . - genotype)
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## unable to evaluate scaled gradient

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
    ##   See ?lme4::convergence and ?lme4::troubleshooting.

``` r
anova(m0, m1, test = "Chisq")
```

    ## Data: cp
    ## Models:
    ## m0: cbind(deaths, at_risk - deaths) ~ tray_size + factor(day) + (1 | tray_id)
    ## m1: cbind(deaths, at_risk - deaths) ~ genotype + tray_size + factor(day) + (1 | tray_id)
    ##    npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
    ## m0   12 283.23 321.55 -129.62    259.23                    
    ## m1   14 284.98 329.68 -128.49    256.98 2.251  2     0.3245

``` r
# Check for overdispersion
rdf <- df.residual(m1)
pearson <- sum(residuals(m1, type = "pearson")^2)
overdispersion_ratio <- pearson / rdf
overdispersion_ratio
```

    ## [1] 1.022025

Check residuals (should be \<2 for good fit)

``` r
m1_Model_DevResid <- resid (m1, type = "deviance")
hist(m1_Model_DevResid)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Check variance of residuals across the range of fitted values (most
points should fall within the grey line)

``` r
m1_Model_Predict <- predict(m1)
binnedplot(m1_Model_Predict, (m1_Model_DevResid))
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Check for outliers (Cook’s distance \>1 indicates outliers)

``` r
plot(m1, which = 4)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Double check fit with DHARMa

``` r
sim <- simulateResiduals(m1)
plot(sim)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
testDispersion(sim)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 1.2062, p-value = 0.352
    ## alternative hypothesis: two.sided

``` r
testUniformity(sim)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  simulationOutput$scaledResiduals
    ## D = 0.075647, p-value = 0.2544
    ## alternative hypothesis: two-sided

``` r
testOutliers(sim)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

    ## 
    ##  DHARMa bootstrapped outlier test
    ## 
    ## data:  sim
    ## outliers at both margin(s) = 0, observations = 180, p-value = 1
    ## alternative hypothesis: two.sided
    ##  percent confidence interval:
    ##  0.000000000 0.008472222
    ## sample estimates:
    ## outlier frequency (expected: 0.00155555555555556 ) 
    ##                                                  0

``` r
plotResiduals(sim, cp$genotype)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plotResiduals(sim, cp$tray_size)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
plotResiduals(sim, cp$day)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

Quantile deviations observed in day, potentially because day is treated
as a factor when there is next to no observable death on two of the
days. Could try a simpler model as random effects have no impact on
variance.

``` r
m1_simpler <- glm(
  cbind(deaths, at_risk - deaths) ~ genotype + tray_size + factor(day),
  data = cp,
  family = binomial(link = "cloglog")
)
summary(m1_simpler)
```

    ## 
    ## Call:
    ## glm(formula = cbind(deaths, at_risk - deaths) ~ genotype + tray_size + 
    ##     factor(day), family = binomial(link = "cloglog"), data = cp)
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -3.2790     0.4083  -8.031 9.64e-16 ***
    ## genotypeCFP-     0.1487     0.3690   0.403 0.686985    
    ## genotypeCFP+    -0.2264     0.3838  -0.590 0.555299    
    ## tray_size50L     0.1974     0.2484   0.795 0.426893    
    ## factor(day)2    -0.7368     0.3358  -2.194 0.028208 *  
    ## factor(day)3    -1.3398     0.4226  -3.170 0.001522 ** 
    ## factor(day)4    -1.6657     0.4855  -3.431 0.000601 ***
    ## factor(day)5   -19.0482  1617.8062  -0.012 0.990606    
    ## factor(day)6    -2.5758     0.7320  -3.519 0.000433 ***
    ## factor(day)7    -1.3149     0.4226  -3.112 0.001859 ** 
    ## factor(day)8    -0.9422     0.3683  -2.558 0.010531 *  
    ## factor(day)9    -1.4381     0.4499  -3.196 0.001392 ** 
    ## factor(day)10  -19.0072  1621.7938  -0.012 0.990649    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 226.23  on 179  degrees of freedom
    ## Residual deviance: 152.81  on 167  degrees of freedom
    ## AIC: 282.98
    ## 
    ## Number of Fisher Scoring iterations: 18

``` r
# Check for overdispersion
rdf <- df.residual(m1_simpler)
pearson <- sum(residuals(m1_simpler, type = "pearson")^2)
overdispersion_ratio <- pearson / rdf
overdispersion_ratio
```

    ## [1] 1.015905

Check model fit with DHARMa

``` r
sim_simpler <- simulateResiduals(m1_simpler)
plot(sim_simpler)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
testDispersion(sim_simpler)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 1.2143, p-value = 0.328
    ## alternative hypothesis: two.sided

``` r
testUniformity(sim_simpler)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  simulationOutput$scaledResiduals
    ## D = 0.064191, p-value = 0.4484
    ## alternative hypothesis: two-sided

``` r
testOutliers(sim_simpler)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

    ## 
    ##  DHARMa bootstrapped outlier test
    ## 
    ## data:  sim_simpler
    ## outliers at both margin(s) = 0, observations = 180, p-value = 1
    ## alternative hypothesis: two.sided
    ##  percent confidence interval:
    ##  0.00000000 0.01111111
    ## sample estimates:
    ## outlier frequency (expected: 0.00216666666666667 ) 
    ##                                                  0

``` r
plotResiduals(sim_simpler, cp$genotype)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plotResiduals(sim_simpler, cp$tray_size)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
plotResiduals(sim_simpler, cp$day)
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

No flagged issues, simpler model seems to fit better.

Calculate hazard ratios based on simple model

``` r
coefs <- summary(m1_simpler)$coefficients
hr <- data.frame(
  term  = rownames(coefs),
  HR    = exp(coefs[, "Estimate"]),
  lower = exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"]),
  upper = exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"]),
  p     = coefs[, "Pr(>|z|)"]
)

hr[grep("genotype|tray_size", hr$term), ]
```

    ##                      term        HR     lower    upper         p
    ## genotypeCFP- genotypeCFP- 1.1603084 0.5629635 2.391479 0.6869847
    ## genotypeCFP+ genotypeCFP+ 0.7974092 0.3758143 1.691956 0.5552991
    ## tray_size50L tray_size50L 1.2182267 0.7485925 1.982489 0.4268929

**Interpretation**

No significant difference in survival between the three genotypes

## Survival predictions at population level, and plotting

As tray size or individual tray had no significant impact in the model,
hazard ratios are pooled for trays in each genotype and day

``` r
# Create a prediction grid for each genotype and day
newd <- expand.grid(
  genotype = levels(cp$genotype),
  tray_size = levels(cp$tray_size),
  day = 1:10
)

newd$day <- factor(newd$day)

# Convert cloglog(eta) to probability of death (hazard) for each day
inv_cloglog <- function(eta) 1 - exp(-exp(eta))

# Predict hazard (population level)
eta_hat <- predict(
  m1_simpler,
  newdata = newd,
  type = "link",
  re.form = NA
)
newd$hazard_hat <- inv_cloglog(eta_hat)

# Average hazards across tray sizes for each genotype × day
haz_by_gen_day <- newd %>%
  group_by(genotype, day) %>%
  summarise(hazard = mean(hazard_hat), .groups = "drop") %>%
  arrange(genotype, as.integer(as.character(day))) %>%
  group_by(genotype) %>%
  mutate(Survival = cumprod(1 - hazard)) %>%
  ungroup() %>%
  mutate(day = as.integer(as.character(day)))
```

Bootstrapping for 95% CIs (x1000)

``` r
stat_fun <- function(fit) {
  eta <- predict(fit, newdata = newd, type = "link", re.form = NA)
  haz <- inv_cloglog(eta)
  tmp <- newd
  tmp$hazard <- haz
  
  out <- tmp %>%
    group_by(genotype, day) %>%
    summarise(hazard = mean(hazard), .groups = "drop") %>%
    arrange(genotype, as.integer(as.character(day))) %>%
    group_by(genotype) %>%
    mutate(Survival = cumprod(1 - hazard)) %>%
    ungroup() %>%
    mutate(day = as.integer(as.character(day)))
  out$Survival
}

set.seed(202)
nsim <- 1000
bt <- bootMer(m1, FUN = stat_fun, nsim = nsim, use.u = FALSE,
              type = "parametric", parallel = "no")

G <- length(levels(cp$genotype))
D <- 10
stopifnot(ncol(bt$t) == G * D)

ci_list <- vector("list", G)
gen_levels <- levels(cp$genotype)

for (g in seq_len(G)) {
  cols <- ((g - 1) * D + 1):(g * D)
  qlo <- apply(bt$t[, cols, drop = FALSE], 2, quantile, probs = 0.05, na.rm = TRUE)
  qhi <- apply(bt$t[, cols, drop = FALSE], 2, quantile, probs = 0.95, na.rm = TRUE)
  ci_list[[g]] <- data.frame(
    genotype = gen_levels[g],
    day = 1:D,
    lo = qlo,
    hi = qhi
  )
}

ci_df <- bind_rows(ci_list)

plot_df <- haz_by_gen_day %>%
  left_join(ci_df, by = c("genotype", "day"))
write.csv(plot_df, "SurvivalPlotdf.csv")
```

Plotting prediction curve with 95% CI ribbon

``` r
genotype_cols2 <- c(
  "CFP+" = "#0B3C8C",
  "CFP-" = "#B11226",
  "G3"   = "#0E7C3A"
)

genotype_cols3 <- c(
  "CFP+" = "#0c4084",
  "CFP-" = "#ac1925",
  "G3"   = "#19a453"
)

genotype_names <- c(
  "CFP+" = expression(vgscKO^"+"),
  "CFP-" = expression(vgscKO^"-"),
  "G3"   = expression(G3)
)

ggplot(plot_df, aes(x = day, y = Survival, color = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = genotype_cols3, labels = genotype_names) +
  scale_fill_manual(values = genotype_cols2, labels = genotype_names) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.8,1)) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Day",
    y = "Predicted survival",
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(color = "Genotype", fill = "Genotype")
```

![](vgscKOLarvalSurvival_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ggsave("LarvalSurvival95CI50btsp.png", width=5, height=5)
```
