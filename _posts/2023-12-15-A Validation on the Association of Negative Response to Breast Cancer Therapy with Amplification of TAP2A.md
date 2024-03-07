---
title: "A Validation on the Association of Negative Response to Breast Cancer Therapy with Amplification of TAP2A"
mathjax: true
layout: post
output: rmarkdown::github_document
categories: project
---
* collaborated with: Weiqiao Shen. Weiqiao contributed to the first half, and my (Ruilin's) work begins from the model fitting section)     

# Introduction
A surrogate endpoint is a clinical trial endpoint used as a substitute for a direct measure of how a patient feels, functions, or survives when the clinical outcomes might take too much time to study or in cases where the clinical benefit of improving the surrogate endpoint is well understood[2]. Validated surrogate endpoints likely provide patients with serious diseases more rapid access to promising therapies. In this project, rather than finding a new one, we decided to validate a previously proved surrogate for breast cancer. 

In 2011, a paper by Desmedt, Christine et al. was published, validating biomarkers predictive of response/resistance to anthracyclines in breast cancer. In that paper, the authors claimed that TOP2A (topoisomerase IIα) amplification was significantly associated with pathological complete response (P ≤ 0.001)[1]. With curiosity, we decided to re-validate this observation on the dataset which the authors published on NCBI (National Center for Biotechnology Information) to answer the question whether TOP2A is truly a predictive biomarker of for breast cancer therapy response.

{% highlight r %}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
{% endhighlight %}


{% highlight r %}
library(mice)
library(missRanger)
library(naniar)
library(tableone)
library(tidyverse)
demo <- read.table("../data/demo.txt", header=TRUE, sep="\t", na.strings = c("NA",""))
demo <- demo %>% 
  mutate(`T` = as.numeric(substr(`T`, 2, 2)),
         N = as.numeric(substr(N, 2, 2))) 

dat <- demo %>% 
  select(-samplename, -filename, -geo_accn) %>%
  filter(FINAL_ANALYSIS == "YES") %>%
  mutate(TOP2Atri = as.factor(TOP2Atri)) %>%
  drop_na(TOP2Atri)

{% endhighlight %}

#### Function Definition Part 1
- Weiqiao contributed to function definition 1 
{% highlight r %}
#random forest imputation
one.rfi <- function(data){
  rf.data <- missRanger(data %>% select(-TOP2Atri, -DMFS_time, -OS_time), pmm.k=3, num.trees=100, verbose=F) %>%
    mutate(TOP2Atri=data$TOP2Atri)
  return(rf.data)
}

#multiple imputation
one.mi <- function(data){
  mi.data <- data %>% 
    select(-TOP2Atri, -DMFS_time, -OS_time) %>% 
    mice(print=F) %>% 
    complete %>% 
    mutate(TOP2Atri=data$TOP2Atri)
  return(mi.data)
}

#extract OR from a glm
#vars: a vector of all predictors used to fit the model
one.glm <- function(data, vars="TOP2Atri"){
  expr <- paste0("pCR ~ ", paste(vars, collapse = " + "))
  coefs <- glm(as.formula(expr), data=data, family = binomial()) %>% coef
  #index <- which(names(coefs) %in% c("TOP2Atri", "TOP2Atri0", "TOP2Atri1"))
  return(coefs %>% exp %>% round(2)) # return rounded OR
}
{% endhighlight %}


# Data
We derived the data from the website of the National Center for Biotechnology Information. The original dataset included clinical information from 120 patients with 166 missing values, among which 28 were exposure (TOP2A amplification status) and 6 were outcome of interest (pathologic complete response status (pCR)). Except for TOP2A and pCR status, the clinical information in this dataset included age, primary tumor size, axillary lymph node status, tumor grade, HER2 status, TOP2A overexpression status, ER status, distant metastasis free survival event and time (days), and overall survival. Among the 120 observations, only 114 of them were defined as evaluable patients and used for final analysis. Among the 114 patients included in the final analysis, 27 had missing values in exposure and none of them had missing values in outcome.

However, after carefully comparing the dataset to Table 1&2 in that paper, we noticed the data they used for final analysis had only 57 missing values (33 in HER2 status and none in exposure or outcome) with a sample size of 139. The researchers applied available case analysis that only focused on the 106 patients after excluding those with missingness in HER2 status. This observation was not consistent with the dataset we retrieved which had missingness in exposure and outcome, and the sample size did not match using available case analysis (see Table 1). Nevertheless, we used what we retrieved from the data source because it had a similar distribution compared to the one that the research team used.

#### Data preparation
{% highlight r %}
da_dat <- demo %>% 
  select(-samplename, -filename, -geo_accn, -FINAL_ANALYSIS) %>% 
  mutate(across(c(agebin, `T`, N, grade, HER2FISHbin, TOP2Atri, ESR1bimod,
                  ERBB2bimod, pCR, DMFS_event, OS_event), as.factor))
CreateTableOne(data=da_dat, strata="TOP2Atri", includeNA=T, test=F)
{% endhighlight %}


# Methods
### Missing Values
We excluded patients with missing values in exposure in data manipulation to avoid imputing biased exposure variables, and applied multiple imputation and random forest imputation separately to fill in missing values. In order to avoid generating direct association between the pCR (response) and TOP2A status (exposure), we removed the exposure from the dataset when imputing missing values. The data imputation involved 13 variables, including age, tumor grade, axillary lymph node status, HER2 status, TOP2A amplification/deletion status, ER status, and indicator of distant metastasis free survival event and overall survival.

#### Missing Value Imputation in R
{% highlight r %}
set.seed(1)
gg_miss_var(dat, show_pct=TRUE)
#miss_var_summary(dat)

#Multiple imputation
#dat.mi <- one.mi(dat)
dat.mi <- one.mi(dat %>% 
                   mutate(topoIHC = ifelse(topoIHC == "<10", 
                                           9, 
                                           topoIHC)) %>% 
                   mutate(topoIHC = as.numeric(topoIHC)))
dat.mi <- within(dat.mi, TOP2Atri <- relevel(TOP2Atri, ref = '0'))
#gg_miss_var(dat.mi, show_pct=TRUE)


#Random forest imputation
dat.rfi <- one.rfi(dat)
dat.rfi <- within(dat.rfi, TOP2Atri <- relevel(TOP2Atri, ref = '0'))

{% endhighlight %}

*Ruilin's contribution begins from here(except for the first chunk under model fitting section `Model Fitting in R`)*

### Model Fitting

Two models were fitted to both imputed datasets generated from multiple imputation and random forest imputation, respectively. Same as the paper, the first model we employed was a simple logistic regression model, which treated pCR as the response variable and TOP2A status as the only covariate. The other model was multiple logistic regression incorporating age (binary) and HER2 (human epidermal growth factor receptor 2) expression level (numeric) in addition to TOP2A. Age was the only available demographic variable, and we thought including age as a covariate could be more helpful compared to only including variables related to the tumor or gene/protein expression levels. HER2 was also included because of its association with breast cancer treatment [1]. Given that TOP2A was coded as a three level factor, odds ratios (OR) for having pCR comparing either 1) TOP2A amplification and normal TOP2A, or 2) TOP2A deletion and normal TOP2A were estimated. Two imputation algorithms were compared through examining log OR estimates generated from these datasets and permutated null distributions.

#### Model Fitting in R
{% highlight r %}

#ORs based on MI data
(actual.OR.mi <- one.glm(dat.mi, vars = c("agebin", "HER2FISH", "TOP2Atri")))
(actual.OR.mi.simple <- one.glm(dat.mi))

#ORs based on RFM data
(actual.OR.rfi <- one.glm(dat.rfi, vars = c("agebin", "HER2FISH", "TOP2Atri")))
(actual.OR.rfi.simple <- one.glm(dat.rfi))

{% endhighlight %}

### Permutation
In order to assess the significance of OR point estimates, permutation tests each with 20000 simulations were performed on 8 estimates (ORs comparing TOP2A levels) generated from 4 models (2 fitted on each dataset generated from each imputation algorithm). The null distribution for ORs were simulated by permuting the response variable, pCR, and followed by refitting the same logistic model (either simple or multiple regression) as they have been fitted before the permutations. 95% CIs and p-values were calculated based on the null distribution. 

#### Function Definitions Part 2
{% highlight r %}

#permutation to assess TOP2A significance using OR
one.perm <- function(data, variables="TOP2Atri") {
  ##permute outcome 
  perm_data <- data
  perm_data$pCR <- data$pCR[sample(nrow(data))]
  ##recalculate ORs
  perm_OR <- one.glm(perm_data, vars  = variables)
  return(perm_OR)
}


#find OR from permutation
perm.OR <- function(data, variables = "TOP2Atri", sim_num = nsim) {
  suppressWarnings(
    replicate(sim_num,
              one.perm(data, variables)))
}


#density plot for permutation with CI and observed value
plot_density <- function(data, var_mi, var_rfi, actual_val,
                         title, ci = c(0.025, 0.975)) {
  ggplot(data) +
    geom_density(aes(x = log(\{\{var_mi\}\}), y = ..density.., color = "mi")) +
    geom_density(aes(x = log(\{\{var_rfi\}\}), y = ..density.., color = "rfi")) +
    labs(title = title,
         x = "log OR") +
    geom_vline(aes(xintercept = log(actual_val), color = "actual")) +
    geom_vline(aes(xintercept = log(quantile(\{\{var_mi\}\}, ci)[1]), color = 'CI')) +
    geom_vline(aes(xintercept = log(quantile(\{\{var_mi\}\}, ci)[2]), color = 'CI')) +
    scale_color_manual(name='Legend',
                       values=c('mi'='darkgreen', 'rfi'='darkblue', 'actual' = 'red', 'CI' = 'orange'))
}


#find CI for the results table
find.CI <- function(data, interval = c(0.025, 0.975)){
  ci <- log(quantile(data, interval)) %>% 
    round(3) %>%
    paste(collapse = ", ")
  return(ci)
}
{% endhighlight %}


# Results
The original paper concluded that TOP2A amplification was significantly associated with pCR. In this project, we reached same conclusions for both the simple,    

```math
\hat{pCR_i} = \hat{\beta_0} + \hat{\beta_1}\text{TOP2A}_i
```

and the multiple logistic regression model,    

```math
\hat{pCR_i} = \hat{\beta_0} + \hat{\beta_1}\text{TOP2A}_i + \hat{\beta_2}\text{age}_i + \hat{\beta_3}\text{HER2}_i
```


where in both models, i represents the $`i^{th}`$ patient in the dataset. The total number of patients being included in this project was 87, who were also included in the data analysis in the original paper and had a valid TOP2A expression level recorded.   

The multiple logistic regression model estimated the log OR of pCR between patients with amplified TOP2A and normal TOP2A being 2.581 (-inf, 1.993). Given that the p-value (0.009) was significant under the 0.05 significance level, we rejected the null hypothesis that the probability of pCR does not associate with TOP2A amplification. Similar conclusion applied to the simple logistic regression, where the log OR was estimated to be 2.725 (-inf, 1.459), with a significant p-value of 0.00015. Table 2 summarized estimates for all four models. 

Models for TOP2A deletion also shared the same conclusion as the original paper: no significance association between TOP2A gene deletion and therapy response. The multiple logistic regression estimated a log OR of 0.824 (-inf, 1.864), with a p-value of 0.191, so we could not reject the null hypothesis that pCR does not associate with TOP2A gene deletion. Similarly, the simple logistic regression had an estimated log OR of 0.116 (-inf, 1.459), with an insignificant p-value (0.116) that was not strong enough to reject the null. 

All 95% CI for log OR under the null distribution shared a common lower bound that approached negative infinity. This was due to the randomness in the permutation process that re-assigned the response variable in a way that the logistic regression model estimated some extreme OR that were infinitely close to 0. The subsequent log transformation expanded the CI to negative infinity. Although the simple and multiple logistic regression models had the same conclusion, estimates from the multiple logistic regression model had wider CI compared to the simpler model. Figure 1 provides visualizations for results in Table 2. 


#### Permutation in R
{% highlight r %}
##Permutation simulations

nsim <- 20000
set.seed(1)

#find adj-OR for two imputed datasets for the model adjusted for age, her2
perm_OR_mi <- perm.OR(dat.mi, c("agebin", "HER2FISH", "TOP2Atri"))
perm_OR_rfi <- perm.OR(dat.rfi, c("agebin", "HER2FISH", "TOP2Atri"))

##put all ORs into a df
all_perm_OR <- as.data.frame(cbind(t(perm_OR_mi), t(perm_OR_rfi)))
names(all_perm_OR)[1:5] <- paste0(names(all_perm_OR)[1:5], "_mi")
names(all_perm_OR)[6:10] <- paste0(names(all_perm_OR)[6:10], "_rfi")


#find crude-OR for the simple logistic regression model
perm_OR_mi_simple <- perm.OR(dat.mi)
perm_OR_rfi_simple <- perm.OR(dat.rfi)
all_perm_OR_simple <- as.data.frame(cbind(t(perm_OR_mi_simple), t(perm_OR_rfi_simple)))
names(all_perm_OR_simple)[1:3] <- paste0(names(all_perm_OR_simple)[1:3], "_mi")
names(all_perm_OR_simple)[4:6] <- paste0(names(all_perm_OR_simple)[4:6], "_rfi")

{% endhighlight %}

#### Dist. after Permutation Visualization in R
{% highlight r %}
##Permutation Visualization
#visualize adjusted model null dist.
plot_density(all_perm_OR, `TOP2Atri-1_mi`, `TOP2Atri-1_rfi`, actual.OR.mi[4],
             "Adjusted log OR between deletion and normal TOP2A from permutation test")

plot_density(all_perm_OR, TOP2Atri1_mi, TOP2Atri1_rfi, actual.OR.mi[5],
             "Adjusted log OR between amplification and normal TOP2A from permutation test")


#visualize unadjusted model null dist.
plot_density(all_perm_OR_simple, `TOP2Atri-1_mi`, `TOP2Atri-1_rfi`, actual.OR.mi.simple[2],
             "Crude log OR between deletion and normal TOP2A from permutation test")

plot_density(all_perm_OR_simple, TOP2Atri1_mi, TOP2Atri1_rfi, actual.OR.mi.simple[3],
             "Crude log OR between amplification and normal TOP2A from permutation test")
{% endhighlight %}


#### Results Compilation in R
{% highlight r %}
#Compile Results into A Table
#all values extracted from the mi dataset for simplicity

#p values
p_val_adj_del <- mean(all_perm_OR$`TOP2Atri-1_mi` > actual.OR.mi[4])
p_val_adj_amp <- mean(all_perm_OR$TOP2Atri1_mi > actual.OR.mi[5])
p_val_crude_del <- mean(all_perm_OR_simple$`TOP2Atri-1_mi` > actual.OR.mi.simple[2])
p_val_crude_amp <- mean(all_perm_OR_simple$TOP2Atri1_mi > actual.OR.mi.simple[3])

#CIs
CI_adj_del <- find.CI(all_perm_OR$`TOP2Atri-1_mi`)
CI_adj_amp <- find.CI(all_perm_OR$TOP2Atri1_mi)
CI_crude_del <- find.CI(all_perm_OR_simple$`TOP2Atri-1_mi`)
CI_crude_amp <- find.CI(all_perm_OR_simple$TOP2Atri1_mi)
  
                          
#compile data into a table
result_table <- data.frame(log(c(actual.OR.mi[4], 
                          actual.OR.mi[5],
                          actual.OR.mi.simple[2],
                          actual.OR.mi.simple[3])),
                        c(p_val_adj_del, 
                          p_val_adj_amp, 
                          p_val_crude_del, 
                          p_val_crude_amp),
                        c(CI_adj_del, 
                          CI_adj_amp, 
                          CI_crude_del, 
                          CI_crude_amp),
                        row.names = c("adjusted model, deletion",
                           "adjusted model, amplification",
                           "crude model, deletion",
                           "crude model, amplification"))
names(result_table) <- c("log OR", "p values", "CI")

knitr::kable(result_table)
{% endhighlight %}

# Discussion 
To conclude, we reached the same results as the original paper had, even though they did not conduct permutation testing. Working beyond their simple logistic regression, we fitted a multiple logistic regression by adjusting age and HER2 expression level. It also gave statistically significant, although slightly more uncertain, conclusions. One possible reason for not getting more precise estimates even though we adjusted for covariates that were known to be correlated with the response could be the variability in imputed values. We did not try complete case analysis to reduce bias. However, imputation itself has a random process, we still introduced additional variation into our datasets. Plus, as the existing values probably were not good representations of the truth, we might exaggerate bias through imputation if complete cases were not representative. Although authors provided descriptive statistics stratifying on response and by all covariates, we unfortunately could not reproduce the whole dataset because the information is limited. So, we could have a conservative conclusion that TOP2A amplification is significantly associated with higher odds of having complete response. To further improve our validation, we need to get the complete dataset and apply the analysis we have done to the complete data. 

In the logistic regression model, we used TOP2A to predict the pathological complete response, adjusted for age and HER2. However, in the original paper, the authors only provided the unadjusted OR. Including additional covariates may also lead to differences in the estimated OR. Additionally, since we used imputation to fill in missing values, we might introduce unwanted association between clinical parameters and outcome of interest. Except for those two concerns, we also need more consideration on the covariates to involve in logistic regression in terms of precision.

### References
1. Desmedt, Christine et al (2011). Multifactorial approach to predicting resistance to 
anthracyclines. Journal of clinical oncology : official journal of the American Society of Clinical Oncology vol. 29,12: 1578-86. doi:10.1200/JCO.2010.31.2231
2. Food and Drug Administration (2017). Surrogate Endpoint Resources for Drug and
Biologic Development. https://www.fda.gov/drugs/development-resources/surrogate-endpoint-resources-drug-and-biologic-development

