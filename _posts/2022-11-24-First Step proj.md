---
title:  "Potential causes for low birthweight"
mathjax: true
layout: post
output: html_document
categories: project
---

{% highlight r %}
knitr::opts_chunk$set(echo = TRUE)
library(rigr)
library(dplyr)
library(Hmisc)
library(knitr)
library(vtable)
{% endhighlight %}

### setup

{% highlight r %}
#read in table
birth_copy <- read.table("BirthsKingCounty2001-Biost514-517-2022.txt", header=T)

#check size 
dim(birth_copy) # dimension: 25 x 17

#converting Y/N character variables to binary variables: sex, smoker, drinker
#if smoker/drinker: 1
#if female: 1
f_converter <- function(x) if_else(x == "N", 0, 1)
birth <- birth_copy %>% 
  mutate_at(c('smoker', 'drinker'), f_converter) %>% 
  mutate(sex = ifelse(sex=='F', 1, 0))

birth$race <- factor(birth$race, levels = c("asian", "black", "hispanic", "white"))
birth$firstep <- factor(birth$firstep)

{% endhighlight %}

# results

univariate summary statistics
{% highlight r %}

descrip(birth, strata = birth$firstep)
#potential confounder: age, married, race (mild),  ?smokeN, ?smoker, education

#histograms
hist.data.frame(birth[,1:5])
hist.data.frame(birth[,6:10])
hist.data.frame(birth[,11:17])

hist.data.frame(birth[birth$firstep==0,1:5])
hist.data.frame(birth[birth$firstep==0,6:10])
hist.data.frame(birth[birth$firstep==0,11:17])

hist.data.frame(birth[birth$firstep==1,1:5])
hist.data.frame(birth[birth$firstep==1,6:10])
hist.data.frame(birth[birth$firstep==1,11:17])

#table 1
#side by side
st(birth, add.median = T, title = "Table 1. Descriptive Statistics of Participants Not Enrolled vs. Enroleed in the FS Program", group = "firstep")

st(birth[birth$firstep==0,], add.median = T, title = "Table 1. Descriptive Statistics of Participants Not Enrolled in the FS Program")

st(birth[birth$firstep==1,], add.median = T, title = "Table 1. Descriptive Statistics of Participants Enrolled in the FS Program")

#histogram for btw
ggplot(birth, aes(x=bwt)) +
  geom_histogram() + 
  ggtitle("Histogram for Birthweight Stratifying on firstep") +
  xlab("Birthweight") + ylab("Count")

#mean difference
mean(birth[birth$firstep==1,]$bwt) # FS
mean(birth[birth$firstep==0,]$bwt) # non-FS
{% endhighlight %}

_outliers_  

The distribution patterns were similar for people who enrolled or not enrolled in the FS program, so for this section we'll look at the total population only.   
By examining table 1 and histograms, smokeN (the number of cigarettes per day during pregnancy) and drinkN (number of drinks per day during pregnancy) had outliers. For smokeN, the mode was at 0 with a 2.633 standard deviation, while outliers were larger than 15. Similar pattern was observed on drinkN, where the mode was at 0, and sd was 0.4756, but 3 out of 2500 mothers drank 10 or more per day. 

_skewness_ 

Birthweights(bwt) showed a slight left skewness. This might be caused by a large number of preterm births. Pre mature babies had smaller weights and the smallest possible weight could be as low as 0g, so the data included many small bwt values. On the other hand, post mature and normal babies may not show obvious overweight, and the heavier half didn't range as broad as the lighter half.   
Other than bwt, gestation (gestational age at birth), wgain (weightgain), parity (prior number of children in family), wpre (mother's prepregnancy weight) and education (years of mother's education) were also highly skewed.  


# Preliminary analysis

{% highlight r %}
t.test(birth[birth$firstep==0,]$bwt, birth[birth$firstep==1,]$bwt)
{% endhighlight %}


Babies whose mothers were not enrolled in the FS program had slightly higher weight than those whose mothers were enrolled.   
By using the formula $\bar{x} \pm 1.96 \times \frac{s}{\sqrt{n}}$,   
the 95% confidence interval for the FS babies is $3359 \pm 1.96 \times \frac{610.8}{\sqrt{403}}$ = $(3299.365, 3418.635)$   
Similarly, the 95% CI for non-FS babies is $(3401.528, 3448.472)$.  
By comparison, we can see that because the non-FS babies had a larger sample size, the confidence interval was narrower than that for FS babies. 

# What variables may cause such a difference in birthweights? 

{% highlight r %}
#t test for weight gain
t.test(birth$wgain[birth$firstep==1], birth$wgain[birth$firstep==0])
{% endhighlight %}

From Table 1, we can conclude that among the 2500 participants, 2097 mothers enrolled in First Steps, and 403 didn't.   
Of those enrolled in FS, 50.62% were single, while only 16.21% of the non-FS mothers were single.   
For the proportion of people on welfare, 5.21% of FS mothers were on welfare, while only 1% of mothers in the non-FS group were on welfare. It's possible that mothers who enrolled in FS tended to need more financial support, which was the result of not having a spouse.   
12.66% of the FS mothers smoked during pregnancy, while only 5.91% of the non-FS mothers smoked.  
FS mothers tend to be (4.67 yrs) younger and (2.37 yrs) less educated than non-FS mothers.  
White mothers (48.7% for FS vs. 71.6% for non-FS) composed the largest proportion for both the FS and non-FS group. A few differences are: the non-FS group had fewest black (5.7%) and Hispanic (6.4%) mothers, and Asian (16.3%) mothers were slightly outnumbered, but were still far fewer than White mothers. The FS group, on the other hand, had a slightly more averaged race distribution. It had 13.5% Asian mothers, 15.3% Black mothers, and 22.4% Hispanic mothers.  (__figure__)  
In terms of amount of maternal weight-gain during pregnancy, no obvious difference was observed between the FS and the non-FS group. The average weight gain for the FS group was 32.988 pounds and 32.143 for the non-FS group. A t-test with p-value = 0.3081 showed that there was no difference among the two. 

### Potential Risk Factors

married, welfare, and smoker are possible dichotomous variables that affect birth weights of the new borns. 
{% highlight r %}
#descriptive stats for the three variables
descrip(birth[c('firstep', 'bwt')], strata = birth$firstep)

##table 2 includes: 1 for FS group, another for non-FS, descriptive data (p-values for t-test, CI, descriptive stats) 

#difference in bwt by firstep
t.test(birth$bwt[birth$firstep==1], birth$bwt[birth$firstep==0])

#t-tests for the three variables
##married
t.test(birth$bwt[birth$married==1], birth$bwt[birth$married==0])

##welfare
t.test(birth$bwt[birth$welfare==1], birth$bwt[birth$welfare==0])

#smoker
t.test(birth$bwt[birth$smoker==1], birth$bwt[birth$smoker==0])

#figure
ggplot(birth, aes(x=married, y=bwt)) +
  geom_boxplot(aes(group=married)) +
  ggtitle("a. Effect of Marital Status on Birthweight") +
  xlab("Marital Status") + ylab("Birthweight") 
ggplot(birth, aes(x=smoker, y=bwt)) +
  geom_boxplot(aes(group=smoker)) +
  ggtitle("b. Effect of Smoking History During Pregnancy on Birthweight") +
  xlab("Smoking During Pregnancy") + ylab("Birthweight") 
ggplot(birth, aes(x=welfare, y=bwt)) +
  geom_boxplot(aes(group=welfare)) +
  ggtitle("c. Effect of Welfare Status on Birthweight") +
  xlab("Welfare Status") + ylab("Birthweight") 
ggplot(birth, aes(x=firstep, y=bwt)) +
  geom_boxplot(aes(group=firstep)) +
  ggtitle("d. Effect of Enrolling in FS on Birthweight") +
  xlab("Enrolled in FS") + ylab("Birthweight") 
{% endhighlight %}

{% highlight r %}
#treatment effect
#diff. between smokers who enrolled in FS vs. not enrolled in FS
diff.smoker = mean(birth$bwt[birth$firstep==1 & birth$smoker==1]) - mean(birth$bwt[birth$firstep==0 & birth$smoker==1])
t.test(birth$bwt[birth$firstep==1 & birth$smoker==1], birth$bwt[birth$firstep==0 & birth$smoker==1])

#diff. between non-smokers who enrolled in FS vs. not enrolled in FS
diff.nonsmoker = mean(birth$bwt[birth$firstep==1 & birth$smoker==0]) - mean(birth$bwt[birth$firstep==0 & birth$smoker==0])
t.test(birth$bwt[birth$firstep==1 & birth$smoker==0], birth$bwt[birth$firstep==0 & birth$smoker==0])
{% endhighlight %}
Among the three risk factors, smoking during pregnancy might be the easiest factor to change. Welfare may get rejected and marital status couldn't be easily changed. 

Smokers who enrolled in FS had babies that were, on average, 67.537g lighter than babies whose mothers didn't enrolled in FS. Non smokers who enrolled in FS had babies that were, on average, 47.99622g lighter than babies whose mothers didn't enrolled in FS. Although non smokers seem to make the FS and non-FS group to have babies with similar weights, both t-tests had p-values larger than 0.05, meaning that the true difference between the two groups are not significantly different.  

{% highlight r %}
#figure stratifying on race
ggplot(birth[!is.na(birth$race),], aes(x=firstep, y=bwt), na.omit(race)) + 
  geom_boxplot() +
  facet_wrap(~race) +
  ggtitle("Boxplots for Birthweights vs. Enrollment in FS, Stratifying on Race") +
  xlab("First Step") + ylab("Birthweight")

ggplot(birth, aes(x=race, y=bwt)) + 
  geom_point() +
  facet_wrap(~firstep)
{% endhighlight %}

  

From the figure above, we may conclude that the birthweights of babies whose mothers enrolled in FS were more centralized around the mean weight. The heaviest FS baby was lighter than the heaviest non-FS baby, and the lightest FS baby was heavier than the non-FS baby in group Asian, Black, and White. For group Hispanic, FS group had a few outliers with extremely low weight babies, but if we exclude those points, the trend in other groups applies to Hispanic group as well, that is, baby weights in FS group had a smaller max weight and larger min weight compared to non-FS group. Strictly speaking, non of the race groups had significant benefit by participating FS, only Hispanic group showed smaller proportion of low weight babies while kept the proportion of normal weight babies when their mothers participated FS. 


# Low Birthweight and Very Low Birthweight

{% highlight r %}
#low BW
low_FS <- nrow(birth[c(birth$firstep==1 & birth$bwt<2500),])/nrow(birth[birth$firstep==1,])
low_nFS <- nrow(birth[c(birth$firstep==0 & birth$bwt<2500),])/nrow(birth[birth$firstep==0,])

#very low BW
vlow_FS <- nrow(birth[c(birth$firstep==1 & birth$bwt<1500),])/nrow(birth[birth$firstep==1,])
vlow_nFS <- nrow(birth[c(birth$firstep==0 & birth$bwt<1500),])/nrow(birth[birth$firstep==0,])

#low BW for var. smoker
low_smoker <- nrow(birth[c(birth$smoker==1 & birth$bwt<2500),])/nrow(birth[birth$smoker==1,])
low_nsmoker <- nrow(birth[c(birth$smoker==0 & birth$bwt<2500),])/nrow(birth[birth$smoker==0,])

#very low BW for smoker
vlow_smoker <- nrow(birth[c(birth$smoker==1 & birth$bwt<1500),])/nrow(birth[birth$smoker==1,])
vlow_nsmoker <- nrow(birth[c(birth$smoker==0 & birth$bwt<1500),])/nrow(birth[birth$smoker==0,])

#low BW for var. married
low_married <- nrow(birth[c(birth$married==1 & birth$bwt<2500),])/nrow(birth[birth$married==1,])
low_nmarried <- nrow(birth[c(birth$married==0 & birth$bwt<2500),])/nrow(birth[birth$married==0,])

#very low BW for var. married
vlow_married <- nrow(birth[c(birth$married==1 & birth$bwt<1500),])/nrow(birth[birth$married==1,])
vlow_nmarried <- nrow(birth[c(birth$married==0 & birth$bwt<1500),])/nrow(birth[birth$married==0,])

#low BW for var. welfare
low_welfare <- nrow(birth[c(birth$welfare==1 & birth$bwt<2500),])/nrow(birth[birth$welfare==1,])
low_nwelfare <- nrow(birth[c(birth$welfare==0 & birth$bwt<2500),])/nrow(birth[birth$welfare==0,])

#very low BW for var. welfare
vlow_welfare <- nrow(birth[c(birth$welfare==1 & birth$bwt<1500),])/nrow(birth[birth$welfare==1,])
vlow_nwelfare <- nrow(birth[c(birth$welfare==0 & birth$bwt<1500),])/nrow(birth[birth$welfare==0,])

{% endhighlight %}

FS group had 6.203% of low BW and 1.489% very low BW babies, while non-FS group only had 4.864% of low BW and 0.525% very low BW babies. FS group had both a higher proportion of low BW babies and very low BW babies, meaning that enrolling in FS didn't help mothers to have stronger babies on average. However, we couldn't conclude that the FS program didn't help to improve babies' weights because we need to consider the group difference between mothers who enrolled in FS and who didn't. It might be that FS mothers would give birth to even weaker babies if they weren't in FS. 

# discussion 

Overall, FS babies did not outperform than non-FS group. 

{% highlight r %}
t.test(birth[birth$married==0 & birth$welfare==1 & birth$smoker==1 & birth$firstep==0,]$bwt, birth[birth$married==0 & birth$welfare==1 & birth$smoker==1 & birth$firstep==1,]$bwt)
       
t.test(birth[birth$married==1 & birth$welfare==0 & birth$smoker==0 & birth$firstep==0,]$bwt, birth[birth$married==1 & birth$welfare==0 & birth$smoker==0 & birth$firstep==1,]$bwt)
{% endhighlight %}

