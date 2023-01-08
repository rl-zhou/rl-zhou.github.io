---
title: "ML model for potential causes of poverty"
mathjax: true
layout: post
output: html_document
categories: project
---

{% highlight r %}
library(tidyverse)
library(fs)
library(dbplyr)
library(maps)
library(stringr)
library(dendextend)
library(rpart.plot)
library(tree)
library(stats)
library(ISLR)
library(class)
library(FNN)
library(purrr)
library(maptree)
library(glmnet)
{% endhighlight %}


{% highlight r %}
state.name <- c(state.name, "District of Columbia")
state.abb <- c(state.abb, "DC")

# read in census data
## exclude Puerto Rico samples to ease visualization 
## transform state names to abbreviations to match the education dataset
census <- read_csv("acs2017_county_data.csv") %>% 
  select(-CountyId, -ChildPoverty, -Income, -IncomeErr, -IncomePerCap, -IncomePerCapErr) %>%
   mutate(State = state.abb[match(`State`, state.name)]) %>%
   filter(State != "PR")
{% endhighlight %}


{% highlight r %}
# read in education data
## change of variables & variable names to match with the census dataset
education <- read_csv("Education.csv") %>%
   filter(!is.na(`2003 Rural-urban Continuum Code`)) %>%
   filter(State != "PR") %>%
   select(-`FIPS Code`,
          -`2003 Rural-urban Continuum Code`,
          -`2003 Urban Influence Code`,
          -`2013 Rural-urban Continuum Code`,
          -`2013 Urban Influence Code`) %>%
   rename(County = `Area name`)
{% endhighlight %}


# Preliminary Data Analysis
### census data

{% highlight r %}
# check that all states are included (should == 51)
length(unique(census$State)) 

# check for missing values
## == 0
sum(is.na(census))
{% endhighlight %}


### education data

{% highlight r %}
# check for missing values
sum(is.na(education))
sum(rowSums(is.na(education))!=0)

# compare the number of distinct counties in two datasets
length(unique(education$County))
length(unique(census$County))
{% endhighlight %}


# Data Wrangling

{% highlight r %}
# exclude missing values, and only include State, County, and the four variables for year 2015-19
edu <- education[complete.cases(education),] %>% 
  select("State","County", "Less than a high school diploma, 2015-19", 
         "High school diploma only, 2015-19","Some college or associate's degree, 2015-19", "Bachelor's degree or higher, 2015-19")  


# total population = the sum of the four columns
edu<- edu%>%
  mutate(total_population = rowSums(edu[, -(1:2)]))


# an aggregated dataset that is grouped based on state
## rename columns
education.state <- edu %>% 
  group_by(State) %>% 
  summarise(`less than high school` = sum(`Less than a high school diploma, 2015-19`),
            `high school` = sum(`High school diploma only, 2015-19`),
            `associate's` = sum(`Some college or associate's degree, 2015-19`), 
            `bachelor's or higher` = sum(`Bachelor's degree or higher, 2015-19`),
            total_population = sum(total_population))


# add a new variable (majority) that records the education level with the largest population in each state
state.level <- education.state %>% 
  mutate(majority = colnames(education.state)[max.col(education.state[,2:5])+1])
{% endhighlight %}


# Visualization

{% highlight r %}
states <- map_data("state")

# merge state.level and states
## not aggregated by states
states <- states %>% 
  mutate(region=state.abb[match(region, tolower(state.name))])
state.merged <- state.level %>% 
  left_join(states, by=c("State" = "region"))


# map the education level (majority)
ggplot(data = state.merged) + 
  geom_polygon(aes(x = long, y = lat, fill = majority, group = group),
               color = "white") + 
  coord_fixed(1.3) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Education Level For the Majority in Each State",
       fill = "Education Level") + 
  scale_fill_discrete(labels=c("Associate's", "Bachelor's or Higher", "High School"))
{% endhighlight %}


### visualization of census data

{% highlight r %}
# merge census with states
## census.state has aggregated population by states
census.state <- census %>% 
  mutate_at(c(6:11, 13:25, 27:31), function(x) x/100*census$TotalPop) %>% 
  group_by(State) %>% 
  summarise(across(TotalPop:Unemployment, list(sum)))

# recalculate the proportions
## find out what's the majority group for each variable category
census.majority <- census.state %>% 
  mutate(gender = colnames(census.state)[max.col(census.state[,3:4])+2]) %>% 
  mutate(race = colnames(census.state)[max.col(census.state[,5:10])+4]) %>% 
  mutate(age = ifelse(VotingAgeCitizen_1>TotalPop_1-VotingAgeCitizen_1, "adult", "non-adults")) %>% 
  mutate(occupation = colnames(census.state)[max.col(census.state[,13:17])+12]) %>%
  mutate(commute = colnames(census.state)[max.col(census.state[,18:22])+17]) %>%
  mutate(work_type = colnames(census.state)[max.col(census.state[,26:29])+25]) %>% 
  select(State, gender:work_type)

# after viewing census.majority, we found that the majority differences only in gender, race and commute(commute only has 1 different value). 
View(census.majority)

# plots
## 
census.long_and_lat <- census.majority %>% 
  left_join(states, by=c("State" = "region")) %>% 
  select(State, gender, race, commute, long:group)

## view by gender
ggplot(data = census.long_and_lat) + 
  geom_polygon(aes(x = long, y = lat, fill = gender, group = group),
               color = "white") + 
  coord_fixed(1.3) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Gender Majority in Each State",
       fill = "Gender") + 
  scale_fill_discrete(labels=c("Men", "Women"))


odds <- function(x) seq_along(x) %% 2 > 0


heatmap(census$Poverty)
{% endhighlight %}



{% highlight r %}
for (i in state.abb){
  a <- census %>%
    filter(State == i)
  b <- ggplot(a, aes(x = Poverty)) +
    geom_histogram(bins = 30, color = "black", fill = "gray") +
    geom_vline(aes(xintercept = mean(Poverty)), 
              linetype = "dashed", size = 0.6) +
    ggtitle("State:", i)
  print(b)
}
{% endhighlight %}


### clean and aggregate census

{% highlight r %}
#excluded NA values
# converted {Men, Employed, VotingAgeCitizen} attributes to percentages
# combined {Hispanic, Black, Native, Asian, Pacific} to {Minority}
# removed {Walk, PublicWork, Construction, Unemployment}; removed {Women}

census.clean <- census[rowSums(is.na(census))==0,] %>% 
  mutate(Men = Men/TotalPop*100, Employed = Employed/TotalPop*100, VotingAgeCitizen = VotingAgeCitizen/TotalPop*100) %>% 
  mutate(Minority = Hispanic+Black+Native+Asian+Pacific) %>% 
  select(-(Hispanic:Pacific), -c(Walk, PublicWork, Construction, Unemployment), -Women)

# check for collinearity: SelfEmployed and WorkAtHome are correlated ()
# mm <- cor(census.clean[,3:ncol(census.clean)])
# which(mm > 0.6, arr.ind = TRUE)
{% endhighlight %}


# Dimensionality reudction: PCA
{% highlight r %}
pca_census <- prcomp(select(census.clean, -c("State", "County")), scale. = TRUE)

# store PC1 and PC2
pc.county <- as.data.frame(pca_census$x[,c("PC1", "PC2")])

# get the loading vectors for PC1 and PC2
pc.rotation <- as.data.frame(pca_census$rotation[,c("PC1", "PC2")])

# get the 3 most important variables in PC1
pc.rotation <- pc.rotation %>% 
  mutate(PC1.abs = abs(PC1), PC2.abs = abs(PC2)) %>% 
  arrange(desc(PC1.abs))
head(pc.rotation, 3)
{% endhighlight %}

  

{% highlight r %}
prvar <- pca_census$sdev^2
pve <- prvar/(sum(prvar))

plot(pve, xlab="Principal Component",
  ylab="Proportion of Variance Explained ", ylim=c(0,1),type='b')

plot(cumsum(pve), xlab="Principal Component ",
  ylab=" Cumulative Proportion of Variance Explained ", ylim=c(0,1), type='b') +
  abline(h=0.9) 

numofPC <-which(cumsum(pve) >= 0.9)[1]
cat("According to the Cumulative Proportion of Variance explained,", numofPC, "PC is needed to capture 90% of the variance for analysis. ")
{% endhighlight %}


# Clusturing
{% highlight r %}
census.dist <- dist(census.clean)
set.seed(1)
census.hclust <- hclust(census.dist, method = "complete")

dendro1 <- as.dendrogram(census.hclust) %>%
  color_branches(., k = 10) %>%
  color_labels(., k = 10) %>%
  set(., "labels_cex", 0.5) %>% 
  set_labels(., labels = census.hclust$labels[order.dendrogram(.)])
  

plot(dendro1, horiz = T, main = "Clustering Using Census.Clean")

pc.county.dist <-dist(pc.county)
set.seed(1)
pc.hclust <- hclust(pc.county.dist, method = "complete")
dendro2 <- as.dendrogram(pc.hclust) %>%
  color_branches(., k = 10) %>%
  color_labels(., k = 10) %>%
  set(., "labels_cex", 0.5) %>%
  set_labels(dendo2, labels = pc.hclust$labels[order.dendrogram(.)])
  
plot(dendro2, horiz = T, main = "Clustering using PC1 and PC2")
{% endhighlight %}


{% highlight r %}
clus = cutree(census.hclust, k=10)
table(clus)
clus2 = cutree(pc.hclust, k=10)
table(clus2)
which(census.clean$County == "Santa Barbara County")
cat("Santa Barbara County resides in cluster ", clus[which(census.clean$County == "Santa Barbara County")], "in
    census.clean constructed hierarchical clustering and in cluster ", clus2[228], " of pc.county clustering. ")

clus[228]

clus1withSB <-census.clean[which(clus == clus[228]), ]
clus2withSB <-census.clean[which(clus2 == clus2[228]),]
{% endhighlight %}

 

{% highlight r %}
all <- census.clean %>%
   left_join(edu, by = c("State"="State", "County"="County")) %>%
   na.omit
{% endhighlight %}


# Modeling
{% highlight r %}
all$Poverty[all$Poverty <= 20] <-0
all$Poverty[!all$Poverty <=20] <-1

set.seed(123)
n <- nrow(all)
idx.tr <- sample.int(n, 0.8*n)
all.tr <- all[idx.tr, ]
all.te <- all[-idx.tr, ]

# all <- all %>%
#   select(-c("TotalPop", "VotingAgeCitizen", "WorkAtHome", "MeanCommute", "Drive", "Transit", "Carpool","OtherTransp", "PrivateWork", "FamilyWork"))
{% endhighlight %}


{% highlight r %}
set.seed(123)
nfold <- 10
folds <- sample(cut(1:nrow(all.tr), breaks=nfold, labels=FALSE))

calc_error_rate = function(predicted.value, true.value){ 
 return(mean(true.value!=predicted.value))
}
records = matrix(NA, nrow=3, ncol=2)
colnames(records) = c("train.error","test.error")
rownames(records) = c("tree","logistic","lasso")
colnames(all.tr)
{% endhighlight %}


{% highlight r %}
colnames(all.tr)[22:25] <- c("LessThanHighSchool", "HighSchool", "Associate", "BachelorHigher")
colnames(all.te)[22:25] <- c("LessThanHighSchool", "HighSchool", "Associate", "BachelorHigher")
{% endhighlight %}


# Classification
### Tree
{% highlight r %}
tree.unpruned = tree(Poverty ~ State + County + TotalPop + Men + Minority + Employed + LessThanHighSchool + HighSchool + Associate + BachelorHigher, data = all.tr)


# tree.unpruned1 = tree(as.factor(Poverty) ~., data = all.tr)


# unpruned predictions
unpruned.tr.pred = predict(tree.unpruned, all.tr, type="class")
unpruned.te.pred = predict(tree.unpruned, all.te, type="class")

# unpruned error rates
## training data
# unpruned.tr.err = table(unpruned.tr.pred, all.tr$Poverty)
# 1-sum(diag(unpruned.tr.err))/sum(unpruned.tr.err)
## test data
# unpruned.te.err = table(unpruned.te.pred, all.te$Poverty)
# 1-sum(diag(unpruned.te.err))/sum(unpruned.te.err)

# plot the unpruned tree
draw.tree(tree.unpruned, nodeinfo=TRUE, cex = 0.4)
title("Poverty Unpruned Tree")

# pruning
cv = cv.tree(tree.unpruned, folds, FUN=prune.misclass, K = nfolds)
best.cv = min(cv$size[cv$dev == min(cv$dev)])
tree.pruned = prune.misclass(tree.unpruned, best=best.cv)

# plot the pruned tree
plot(tree.pruned)
text(tree.pruned, pretty=0, col = "blue", cex = .5)
title("Poverty pruned tree of size 3")

# pruned predictions
pruned.tr.pred = predict(tree.pruned, newdata = all.tr, type="class")
pruned.te.pred = predict(tree.pruned, newdata = all.te, type="class")

# pruned error rates
## training data
pruned.tr.err = table(pruned.tr.pred, all.tr$Poverty)
records[1,1] = 1-sum(diag(pruned.tr.err))/sum(pruned.tr.err)
## test data
pruned.te.err = table(pruned.te.pred, all.te$Poverty)
records[1,2] = 1-sum(diag(pruned.te.err))/sum(pruned.te.err)
records
{% endhighlight %}
   


### Logistic

{% highlight r %}
# model
glm.fit = glm(Poverty~., data = all.tr, family=binomial)
summary(glm.fit)

prob.tr = predict(glm.fit, type="response")
logistic.pred.tr=as.factor(ifelse(prob.tr >0.5, 1, 0))
glm.tr.table = table(pred=logistic.pred.tr, true=all.tr$Poverty) # training error table


prob.te = predict(glm.fit, type="response", newdata = all.te)
logistic.pred.te=as.factor(ifelse(prob.te >0.5, 1, 0))
glm.te.table = table(pred=logistic.pred.te, true=all.te$Poverty) # test error table
records[2,1] = 1-sum(diag(glm.tr.table)/sum(glm.tr.table))
records[2,2] = 1-sum(diag(glm.te.table)/sum(glm.te.table))
{% endhighlight %}


### Lasso

{% highlight r %}
set.seed(123)
x = model.matrix(Poverty~., all.tr)[,-1]
x.test = model.matrix(Poverty~., all.te)[,-1]
y = all.tr$Poverty
y.test = all.te$Poverty

lasso.mod<-glmnet(x, y, alpha = 1, family = 'binomial', lambda = seq(1,20)*1e-5)
cv.out.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial")
bestlam = cv.out.lasso$lambda.min
lasso.pred.test = predict(lasso.mod, s= bestlam, newx = x.test,
                          type ="response")
lasso.pred.train = predict(lasso.mod, s= bestlam, newx = x,
                          type ="response")

lasso.pred.train <- ifelse(lasso.pred.train > .5, 1, 0)
lasso.pred.test <- ifelse(lasso.pred.test > .5, 1, 0)

lasso.coef = predict(lasso.mod, type = "coefficients", s = bestlam)


latrain <- table(pred= lasso.pred.train, true = y)
latest <- table(pred= lasso.pred.test, true = y.test)
latest
latrain

plot(lasso.mod, xvar="lambda", label = TRUE)
plot(cv.out.lasso)
abline(v = log(cv.out.lasso$lambda.min), col="red", lwd=3, lty=2)

records[3,1] = 1-sum(diag(latrain)/sum(latrain))
records[3,2] = 1-sum(diag(latest)/sum(latest))
{% endhighlight %}


{% highlight r %}
records
{% endhighlight %}


### ROC 

{% highlight r %}
# tree
library(ROCR)
tr.pred <-predict(tree.pruned, newdata = all.te)
pred.tree = prediction(tr.pred[,2], all.te$Poverty)
perf.tree = performance(pred.tree, measure="tpr", x.measure="fpr")

# logistic
pred.logistic = prediction(prob.tr, all.tr$Poverty)
perf.logistic = performance(pred.logistic, measure="tpr", x.measure="fpr")


# lasso
pred.lasso = prediction(lasso.pred.test, y.test)
perf.lasso = performance(pred.logistic, measure = "tpr", x.measure="fpr")

plot(perf.tree, lwd = 0.6, col = 2)
plot(perf.logistic,  lwd = 0.6, add = TRUE, col =3)
plot(perf.lasso,  lwd = 0.6,  add = TRUE,col = 4)
legend(0.6,0.6, c("Decision Tree", "Logistic Classification", "Lasso classification"), 2:4)
{% endhighlight %}



### KNN
{% highlight r %}

YTrain = all.tr$Poverty

XTrain = all.tr %>% 
  select(-Poverty) %>%
  scale(., center = TRUE, scale = TRUE)
# XTrain
YTest = all.te$Poverty

XTest = all.te %>% 
  select(-Poverty) %>%
  scale(., center = TRUE, scale = TRUE)

set.seed(123)
cl=as.vector(all.tr$Poverty)
pred.TrainY = knn(train = XTrain, test = XTrain, cl, k=10)

conf.train = table(predicted=pred.TrainY, true = YTrain)
conf.train
TrainerrorKNN <- 1 - sum(diag(conf.train)/sum(conf.train))
TrainerrorKNN

cl2 = as.vector(all.te$Poverty)
pred.TestY = knn(train = XTest, test = XTest, cl = cl2, k=10)
conf.test = table(predicted=pred.TestY, true = YTest)
TesterrorKNN <-1-sum(diag(conf.test)/sum(conf.test))
TesterrorKNN
cat("The trainning error rate for KNN is ", TrainerrorKNN, ", and the test error rate is ", TesterrorKNN, ".")
{% endhighlight %}

### Random Forest
{% highlight r %}
set.seed(123)
# Installing package
# install.packages("caTools")       # For sampling the dataset
# install.packages("randomForest")  # For implementing random forest algorithm

# Loading package
library(caTools)
library(randomForest)
  
# Splitting data in train and test data
# train and test data using the same as the KNN
RFmodel <- randomForest(x=XTrain, 
                        y=as.factor(YTrain), ntree=500, importance =T)
RFmodel
RF_ytr = predict(RFmodel, newdata = XTrain)
# compute the error rates
conf.RF.train = table(predict=RF_ytr, true = YTrain)
conf.trainRF = 1-sum(diag(conf.RF.train)/sum(conf.RF.train))
conf.trainRF
RF_ypred = predict(RFmodel, newdata = XTest)
# compute the error rates
conf.RF.test = table(predict=RF_ypred, true = YTest)
conf.testRF = 1-sum(diag(conf.RF.test)/sum(conf.RF.test))
conf.testRF

cat("The trainning error rate for randomforest is ", conf.trainRF, " and the test error rate is ", conf.testRF, ".")
varImpPlot(RFmodel, sort=T, main="Variable Importance for RF-poverty vs others", n.var=5)
# list out the 5 most important variables recognized by the RF through constructing 500 trees
{% endhighlight %}

 

{% highlight r %}
Poverty.num <- census %>%
   left_join(edu, by = c("State"="State", "County"="County")) %>%
  na.omit

colnames(Poverty.num)[32:36] <- c("LessThanHighSchool", "HighSchool", "Associate", "BachelorHigher")
colnames(Poverty.num)[32:36] <- c("LessThanHighSchool", "HighSchool", "Associate", "BachelorHigher")

set.seed(123)
nl <- nrow(all)
idx.trl <- sample.int(nl, 0.8*nl)
Poverty.num.tr <- Poverty.num[idx.tr, ]
Poverty.num.te <- Poverty.num[-idx.tr, ]

# lm.mod <- lm(Poverty ~.-(County + State), data = Poverty.num.tr)
# summary(lm.mod)
# select significant variables
lm.mod <- lm(Poverty ~ Black + Native + VotingAgeCitizen + MeanCommute + Employed + Unemployment + HighSchool + Associate , data=Poverty.num.tr)
# prediction raw value
lm.tr.pred <- predict(lm.mod, newdata = Poverty.num.tr[,-14], type = "response")
lm.te.pred <- predict(lm.mod, newdata = Poverty.num.te[,-14], type = "response")

# convert prediction back to 0&1
lm.tr.pred <- ifelse(lm.tr.pred<=20, 0, 1)
lm.te.pred <- ifelse(lm.te.pred<=20, 0, 1)

# error rate
lm.tr.err.table <- table(pred=lm.tr.pred, true = all.tr$Poverty)
lm.te.err.table <- table(pred=lm.te.pred, true = all.te$Poverty)

tr.err.rate <- 1-sum(diag(lm.tr.err.table)/sum(lm.tr.err.table))
te.err.rate <- 1-sum(diag(lm.te.err.table)/sum(lm.te.err.table))

tr.err.rate
te.err.rate
{% endhighlight %}



{% highlight r %}
summary(lm.mod)
plot(lm.mod$residuals ~ lm.mod$fitted.values, xlab = "Fitted", ylab = "Residual")
title('Residual vs. fitted plot')
plot(predict(lm.mod), Poverty.num.tr$Poverty, xlab = "Predicted Values", ylab = "Observed Values")
title("Observed vs. Predicted")
abline(a = 0, b = 1, col = "red", lwd = 2)

# tree
tr.pred <-predict(tree.pruned, newdata = all.te)
pred.tree = prediction(tr.pred[,2], all.te$Poverty)
perf.tree = performance(pred.tree, measure="tpr", x.measure="fpr")

# logistic
pred.logistic = prediction(prob.tr, all.tr$Poverty)
perf.logistic = performance(pred.logistic, measure="tpr", x.measure="fpr")

# lasso
pred.lasso = prediction(lasso.pred.test, y.test)
perf.lasso = performance(pred.logistic, measure = "tpr", x.measure="fpr")

#linear model
pred.lm = prediction(lm.te.pred, all.te$Poverty)
perf.lm = performance(pred.lm, measure = "tpr", x.measure ="fpr")



plot(perf.tree, lwd = 0.6, col = 2)
plot(perf.logistic,  lwd = 0.6, add = TRUE, col =3)
plot(perf.lasso,  lwd = 0.6,  add = TRUE,col = 4)
plot(perf.lm, lwd = 0.6, add = TRUE, col =5)

legend(0.6,0.6, c("Decision Tree", "Logistic Classification", "Lasso classification", "Linear Regression->classification"), 2:5)
{% endhighlight %}



