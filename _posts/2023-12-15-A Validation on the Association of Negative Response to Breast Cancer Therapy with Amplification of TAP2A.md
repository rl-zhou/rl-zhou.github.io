---
title: "A Validation on the Association of Negative Response to Breast Cancer Therapy with Amplification of TAP2A"
mathjax: true
layout: post
output: html_document
categories: project
---
* collaborated with: Weiqiao Shen

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

### random forest imputation
{% highlight r %}
one.rfi <- function(data){
  rf.data <- missRanger(data %>% select(-TOP2Atri, -DMFS_time, -OS_time), pmm.k=3, num.trees=100, verbose=F) %>%
    mutate(TOP2Atri=data$TOP2Atri)
  return(rf.data)
}

# multiple imputation
one.mi <- function(data){
  mi.data <- data %>% 
    select(-TOP2Atri, -DMFS_time, -OS_time) %>% 
    mice(print=F) %>% 
    complete %>% 
    mutate(TOP2Atri=data$TOP2Atri)
  return(mi.data)
}

# extract OR from a glm
# vars: a vector of all predictors used to fit the model
one.glm <- function(data, vars="TOP2Atri"){
  expr <- paste0("pCR ~ ", paste(vars, collapse = " + "))
  coefs <- glm(as.formula(expr), data=data, family = binomial()) %>% coef
  #index <- which(names(coefs) %in% c("TOP2Atri", "TOP2Atri0", "TOP2Atri1"))
  return(coefs %>% exp %>% round(2)) # return rounded OR
}


# permutation to assess TOP2A significance using OR
one.perm <- function(data, variables="TOP2Atri") {
  # permute outcome 
  perm_data <- data
  perm_data$pCR <- data$pCR[sample(nrow(data))]
  
  # recalculate ORs
  perm_OR <- one.glm(perm_data, vars  = variables)
  
  return(perm_OR)
}

# bootstrap
# impute each bootstrap
one.boot <- function(data, method, variables = "TOP2Atri"){
  
  # bootstrap data
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  
  # imputation
  if (method == "rfi") {
    tmp.dat <- one.rfi(boot_data)
  } else if (method == "mi") {
    tmp.dat <- one.mi(boot_data %>% 
                   mutate(topoIHC = ifelse(topoIHC == "<10", 
                                           9, 
                                           topoIHC)) %>% 
                   mutate(topoIHC = as.numeric(topoIHC)))
  }
  
  tmp.dat <- within(tmp.dat, TOP2Atri <- relevel(TOP2Atri, ref = '0'))
  
  # resampled OR
  boot_OR <- one.glm(boot_data, vars = variables)
  
  return(boot_OR)
}

tmp.one.boot <- function(data, variables = "TOP2Atri"){
  
  # bootstrap data
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  
  # resampled OR
  boot_OR <- one.glm(boot_data, vars = variables)
  
  return(boot_OR)
}


# find log CI
get.CI <- function(sampling_dist, CI_range = c(0.025,0.975), actual_OR){
  quants <- quantile(sampling_dist, CI_range)
  # dists <- quants - actual_OR
  # CI <- actual_OR - dists[2:1]
  upper_dist <- quants[2] - actual_OR# the upper distance
  lower_dist <- actual_OR - quants[1] # the lower distance
  CI <- c(actual_OR - upper_dist, actual_OR + lower_dist) # pivoting the confidence interval
  names(CI) <- names(quants)

  if (CI[1] < 0) {
    CI[1] = 0
  }
  return(log(CI))
}

{% endhighlight %}

{% highlight r %}
da_dat <- demo %>% 
  select(-samplename, -filename, -geo_accn, -FINAL_ANALYSIS) %>% 
  mutate(across(c(agebin, `T`, N, grade, HER2FISHbin, TOP2Atri, ESR1bimod,
                  ERBB2bimod, pCR, DMFS_event, OS_event), as.factor))
CreateTableOne(data=da_dat, strata="TOP2Atri", includeNA=T, test=F)
{% endhighlight %}



{% highlight r %}
# Should exposure be used to impute outcome? If yes, do we need to apply cross validation in imputation?
# How to impute exposure variable? or simply removing missing values in exposure?

set.seed(1)
gg_miss_var(dat, show_pct=TRUE)
#miss_var_summary(dat)

# Multiple imputation
#dat.mi <- one.mi(dat)
dat.mi <- one.mi(dat %>% 
                   mutate(topoIHC = ifelse(topoIHC == "<10", 
                                           9, 
                                           topoIHC)) %>% 
                   mutate(topoIHC = as.numeric(topoIHC)))
dat.mi <- within(dat.mi, TOP2Atri <- relevel(TOP2Atri, ref = '0'))
# gg_miss_var(dat.mi, show_pct=TRUE)


# Random forest imputation
dat.rfi <- one.rfi(dat)
dat.rfi <- within(dat.rfi, TOP2Atri <- relevel(TOP2Atri, ref = '0'))

# Data based on table 2
# Amplification and non-amplification only
dat.tbl2.1 <- data.frame(pCR=c(rep(0, 96), rep(1, 10)),
                        TOP2Atri=c(rep(1, 9), 
                                   rep(0, 96-9),
                                   rep(1, 6),
                                   rep(0, 10-6)))

{% endhighlight %}


{% highlight r %}

# ORs based on MI data
(actual.OR.mi <- one.glm(dat.mi, vars = c("agebin", "HER2FISH", "TOP2Atri")))
(actual.OR.mi.simple <- one.glm(dat.mi))

# ORs based on RFM data
(actual.OR.rfi <- one.glm(dat.rfi, vars = c("agebin", "HER2FISH", "TOP2Atri")))
(actual.OR.rfi.simple <- one.glm(dat.rfi))

# ORs based on dat.tbl2.1
(actual.OR.tbl2 <- one.glm(dat.tbl2.1))

# ORs based on dat.tbl2.2
#(actual.OR.tbl2 <- one.glm(dat.tbl2.2))

{% endhighlight %}




{% highlight r %}

# run simulations
nsim <- 20000
set.seed(1)

# adjust for age, her2
perm_OR_mi <- suppressWarnings(replicate(nsim, one.perm(dat.mi, variables = c("agebin", "HER2FISH", "TOP2Atri"))))
perm_OR_rfi <- suppressWarnings(replicate(nsim, one.perm(dat.rfi, variables = c("agebin", "HER2FISH", "TOP2Atri"))))
all_perm_OR <- as.data.frame(cbind(t(perm_OR_mi), t(perm_OR_rfi)))
names(all_perm_OR)[1:5] <- paste0(names(all_perm_OR)[1:5], "_mi")
names(all_perm_OR)[6:10] <- paste0(names(all_perm_OR)[6:10], "_rfi")

# not adjust anything: simple logistic
perm_OR_mi_simple <- suppressWarnings(replicate(nsim, one.perm(dat.mi)))
perm_OR_rfi_simple <- suppressWarnings(replicate(nsim, one.perm(dat.rfi)))
all_perm_OR_simple <- as.data.frame(cbind(t(perm_OR_mi_simple), t(perm_OR_rfi_simple)))
names(all_perm_OR_simple)[1:3] <- paste0(names(all_perm_OR_simple)[1:3], "_mi")
names(all_perm_OR_simple)[4:6] <- paste0(names(all_perm_OR_simple)[4:6], "_rfi")

# table 2
perm_OR_tbl2 <- suppressWarnings(replicate(nsim, one.perm(dat.tbl2.1)))
perm_OR_tbl2 <- as.data.frame(t(perm_OR_tbl2))


# adjusted model
ggplot() + 
  geom_density(data = all_perm_OR, aes(x = log(`TOP2Atri-1_mi`), y = ..density.., color = "mi")) + 
  geom_density(data = all_perm_OR, aes(x = log(`TOP2Atri-1_rfi`), y = ..density.., color = "rfi")) + 
  labs(x = "log OR",
       title = "Adjusted log OR between deletion and normal TOP2A from permutation test") +
  geom_vline(aes(xintercept = log(actual.OR.mi[4]), color = "actual")) + 
  geom_vline(aes(xintercept = log(quantile(all_perm_OR$`TOP2Atri-1_mi`, 
                                       c(0.025, 0.975))),
                 color = 'CI')) + 
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'actual' = 'red', 'CI' = 'orange'))

ggplot() + 
  geom_density(data = all_perm_OR, aes(x = log(TOP2Atri1_mi), y = ..density.., color = "mi")) + 
  geom_density(data = all_perm_OR, aes(x = log(TOP2Atri1_rfi), y = ..density.., color = "rfi")) + 
  labs(title = "Adjusted log OR between amplification and normal TOP2A from permutation test",
       x = "log OR") +
  geom_vline(aes(xintercept = log(actual.OR.mi[5]), color = "actual")) +
  geom_vline(aes(xintercept = log(quantile(all_perm_OR$TOP2Atri1_rfi, 
                                       c(0.025, 0.975))),
                 color = 'CI')) +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'actual' = 'red', 'CI' = 'orange'))


# unadjusted model
ggplot() + 
  geom_density(data = all_perm_OR_simple, aes(x = log(`TOP2Atri-1_mi`), y = ..density.., color = "mi")) + 
  geom_density(data = all_perm_OR_simple, aes(x = log(`TOP2Atri-1_rfi`), y = ..density.., color = "rfi")) + 
  labs(x = "log OR",
       title = "Crude log OR between deletion and normal TOP2A from permutation test") +
  geom_vline(aes(xintercept = log(actual.OR.mi.simple[2]), color = "actual")) +
  geom_vline(aes(xintercept = log(quantile(all_perm_OR$`TOP2Atri-1_rfi`, 
                                       c(0.025, 0.975))),
                 color = 'CI')) +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'actual' = 'red', 'CI' = 'orange'))

ggplot() + 
  geom_density(data = all_perm_OR_simple, aes(x = log(TOP2Atri1_mi), y = ..density.., color = "mi")) + 
  geom_density(data = all_perm_OR_simple, aes(x = log(TOP2Atri1_rfi), y = ..density.., color = "rfi")) + 
  labs(title = "Crude log OR between amplification and normal TOP2A from permutation test",
       x = "log OR") +
  geom_vline(aes(xintercept = log(actual.OR.mi.simple[3]), color = "actual")) +
  geom_vline(aes(xintercept = log(quantile(all_perm_OR$TOP2Atri1_rfi, 
                                       c(0.025, 0.975))),
                 color = 'CI')) +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'actual' = 'red', 'CI' = 'orange'))

# table 2
ggplot() + 
  geom_density(data = perm_OR_tbl2, aes(x = log(TOP2Atri), y = ..density..), color = "darkblue") + 
  labs(x = "Crude log OR between amplification and non-amplification TOP2A") +
  geom_vline(aes(xintercept = log(actual.OR.mi.simple[2])), color = "red") +
  geom_vline(aes(xintercept = log(quantile(perm_OR_tbl2$TOP2Atri, 
                                       c(0.025, 0.975))),
                 color = 'miCI')) 

{% endhighlight %}




{% highlight r %}

set.seed(1)
nsim <- 1000

boot_OR_mi <- suppressWarnings(replicate(nsim, one.boot(dat, "mi", variables = c("agebin", "HER2FISH", "TOP2Atri"))))
boot_OR_rfi <- suppressWarnings(replicate(nsim, one.boot(dat, "rfi", variables = c("agebin", "HER2FISH", "TOP2Atri"))))
all_boot_OR <- as.data.frame(cbind(t(boot_OR_mi), t(boot_OR_rfi)))
names(all_boot_OR) <- names(all_perm_OR)

boot_OR_mi_simple <- suppressWarnings(replicate(nsim, one.boot(dat, "mi")))
boot_OR_rfi_simple <- suppressWarnings(replicate(nsim, one.boot(dat, "rfi")))
all_boot_OR_simple <- as.data.frame(cbind(t(boot_OR_mi_simple), t(boot_OR_rfi_simple)))
names(all_boot_OR_simple) <- names(all_perm_OR_simple)


# visualization
## adjusted model
ggplot() + 
  geom_density(data = all_boot_OR, aes(x = log(`TOP2Atri-1_mi`), y = ..density.., color = "mi")) + 
  geom_density(data = all_boot_OR, aes(x = log(`TOP2Atri-1_rfi`), y = ..density.., color = "rfi")) + 
  labs(x = "Crude log OR between deletion and normal TOP2A",
       title = "Bootstrapping Adjusted log OR with Imputations, Deletion") +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR$`TOP2Atri-1_mi`, 
                                      actual_OR = actual.OR.mi.simple[2])), 
                 color = 'miCI')) +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR$`TOP2Atri-1_rfi`, 
                                      actual_OR = actual.OR.rfi.simple[2])), color = 'rfiCI')) +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'miCI' = 'orange', 'rfiCI' = 'darkred'))

ggplot() + 
  geom_density(data = all_boot_OR, aes(x = log(TOP2Atri1_mi), y = ..density.., color = "mi")) + 
  geom_density(data = all_boot_OR, aes(x = log(TOP2Atri1_rfi), y = ..density.., color = "rfi")) + 
  geom_vline(aes(xintercept = (get.CI(all_boot_OR$TOP2Atri1_mi, 
                                      actual_OR = actual.OR.mi.simple[3])), 
                 color = 'miCI')) +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR$TOP2Atri1_rfi, 
                                      actual_OR = actual.OR.rfi.simple[3])), 
                 color = 'rfiCI')) +
  labs(x = "Crude log OR between amplification and normal TOP2A",
       title = "Bootstrapping Adjusted log OR with Imputations, Amplification") +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'miCI' = 'orange', 'rfiCI' = 'darkred'))


## crude model
ggplot() + 
  geom_density(data = all_boot_OR_simple, aes(x = log(`TOP2Atri-1_mi`), y = ..density.., color = "mi")) + 
  geom_density(data = all_boot_OR_simple, aes(x = log(`TOP2Atri-1_rfi`), y = ..density.., color = "rfi")) + 
  labs(x = "Crude log OR between deletion and normal TOP2A",
       title = "Bootstrapping Crude log OR with Imputations, Deletion") +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR_simple$`TOP2Atri-1_mi`, 
                                      actual_OR = actual.OR.mi.simple[2])), 
                 color = 'miCI')) +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR_simple$`TOP2Atri-1_rfi`, 
                                      actual_OR = actual.OR.rfi.simple[2])), color = 'rfiCI')) +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'miCI' = 'orange', 'rfiCI' = 'darkred'))

ggplot() + 
  geom_density(data = all_boot_OR_simple, aes(x = log(TOP2Atri1_mi), y = ..density.., color = "mi")) + 
  geom_density(data = all_boot_OR_simple, aes(x = log(TOP2Atri1_rfi), y = ..density.., color = "rfi")) + 
  geom_vline(aes(xintercept = (get.CI(all_boot_OR_simple$TOP2Atri1_mi, 
                                      actual_OR = actual.OR.mi.simple[3])), 
                 color = 'miCI')) +
  geom_vline(aes(xintercept = (get.CI(all_boot_OR_simple$TOP2Atri1_rfi, 
                                      actual_OR = actual.OR.rfi.simple[3])), 
                 color = 'rfiCI')) +
  labs(x = "Crude log OR between amplification and normal TOP2A",
       title = "Bootstrapping Crude log OR with Imputations, Amplification") +
  scale_color_manual(name='Legend',
                     values=c('mi'='darkgreen', 'rfi'='darkblue', 'miCI' = 'orange', 'rfiCI' = 'darkred'))
{% endhighlight %}





{% highlight r %}
# all values extracted from the mi dataset for simplicity

# p values
p_val_adj_del <- mean(all_perm_OR$`TOP2Atri-1_mi` > actual.OR.mi[4])
p_val_adj_amp <- mean(all_perm_OR$TOP2Atri1_mi > actual.OR.mi[5])
p_val_crude_del <- mean(all_perm_OR_simple$`TOP2Atri-1_mi` > actual.OR.mi.simple[2])
p_val_crude_amp <- mean(all_perm_OR_simple$TOP2Atri1_mi > actual.OR.mi.simple[3])

# CIs
CI_adj_del <- log(quantile(all_perm_OR$`TOP2Atri-1_mi`, 
                                       c(0.025, 0.975))) %>% 
                            round(3) %>% 
                            paste(collapse = ", ")
CI_adj_amp <- log(quantile(all_perm_OR$TOP2Atri1_mi, 
                                       c(0.025, 0.975))) %>% 
                            round(3) %>% 
                            paste(collapse = ", ")
CI_crude_del <- log(quantile(all_perm_OR_simple$`TOP2Atri-1_mi`, 
                                       c(0.025, 0.975))) %>% 
                            round(3) %>% 
                            paste(collapse = ", ")

CI_crude_amp <- log(quantile(all_perm_OR_simple$TOP2Atri1_mi,
                     c(0.025, 0.975))) %>% 
                            round(3) %>% 
                            paste(collapse = ", ")
  
                          

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


