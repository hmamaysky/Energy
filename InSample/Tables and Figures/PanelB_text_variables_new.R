setwd("C:/Users/j1kxb09/Downloads")

# Libraries 
library(sjPlot)
library(lmtest)
library(dplyr)
library(tidyr)
library(sandwich)
library(readstata13)
library(selectiveInference)
library(imputeTS)

data.prices <- read.dta13("transformed_data_prices_v19.2_mod.dta") 
#data.physical <- read.dta13("transformed_data_physical_v19.2_mod.dta")

##########PanelB#############
x <- data.prices[,"PCAsent_Fri"]
x <- na.omit(x)
mean(x)
sd(x)
quantile(x, probs = c(.05, .5, .95))


c <- data.prices[,"PCAfreq_Fri"]
c <- na.omit(c)
mean(c)
sd(c)
quantile(c, probs = c(.05, .5, .95))


d <- data.prices[,"PCAall_Fri"]
d <- na.omit(d)
mean(d)
sd(d)
quantile(d, probs = c(.05, .5, .95))


e <- data.prices[,"artcount_4wk_Fri"]
e <- na.omit(e)
mean(e)
sd(e)
quantile(e, probs = c(.05, .5, .95))

#entropy
f <- data.prices[,"entropy_4wk_Fri"]
f <- na.omit(f)
mean(f)
sd(f)
quantile(f, probs = c(.05, .5, .95))

#######OtherVariables####

#sCo
g <- data.prices[,"stopic1_4wk_Fri"] 
g <- na.omit(g)
mean(g)
sd(g)
quantile(g, probs = c(.05, .5, .95))

#fCo
h <- data.prices[,"ftopic1_4wk_Fri"]
h <- na.omit(h)
mean(h)
sd(h)
quantile(h, probs = c(.05, .5, .95))

#sGom
i <- data.prices[,"stopic2_4wk_Fri"]
i <- na.omit(i)
mean(i)
sd(i)
quantile(i, probs = c(.05, .5, .95))

#fGom
j <- data.prices[,"ftopic2_4wk_Fri"]
j <- na.omit(j)
mean(j)
sd(j)
quantile(j, probs = c(.05, .5, .95))

#sEnv
k <- data.prices[,"stopic3_4wk_Fri"]
k <- na.omit(k)
mean(k)
sd(k)
quantile(k, probs = c(.05, .5, .95))

#fEnv
l <- data.prices[,"ftopic3_4wk_Fri"]
l <- na.omit(l)
mean(l)
sd(l)
quantile(l, probs = c(.05, .5, .95))


#sEpg
m <- data.prices[,"stopic4_4wk_Fri"]
m <- na.omit(m)
mean(m)
sd(m)
quantile(m, probs = c(.05, .5, .95))

#fEpg
n <- data.prices[,"ftopic4_4wk_Fri"]
n <- na.omit(n)
mean(n)
sd(n)
quantile(n, probs = c(.05, .5, .95))


#sBbl
o <- data.prices[,"stopic5_4wk_Fri"]
o <- na.omit(o)
mean(o)
sd(o)
quantile(o, probs = c(.05, .5, .95))


#fBbl
p <- data.prices[,"ftopic5_4wk_Fri"]
p <- na.omit(p)
mean(p)
sd(p)
quantile(p, probs = c(.05, .5, .95))


#sRpc
q <- data.prices[,"stopic6_4wk_Fri"]
q <- na.omit(q)
mean(q)
sd(q)
quantile(q, probs = c(.05, .5, .95))


#fRpc
r <- data.prices[,"ftopic6_4wk_Fri"]
r <- na.omit(r)
mean(r)
sd(r)
quantile(r, probs = c(.05, .5, .95))


#sEp
s <- data.prices[,"stopic7_4wk_Fri"]
s <- na.omit(s)
mean(s)
sd(s)
quantile(s, probs = c(.05, .5, .95))


#fEp
t <- data.prices[,"ftopic7_4wk_Fri"]
t <- na.omit(t)
mean(t)
sd(t)
quantile(t, probs = c(.05, .5, .95))




