library(data.table) #converts data to time series format
library(ggplot2)#to plot various plots
library(fpp2)#examine seasonality graphically
library(forecast) #for various functions related to time series
library(stats)#for acf, Ljung-Box tests
library(tseries)#for applying Dickey Fuller tests for checking stationary
library(here)
library(tidyverse)
library(tsutils)
library(TSA)
library(zoo)
library(TSstudio)
library(Rssa)
library(astsa)
library(urca)
library(quantmod)
library(vrtest)
library(dplyr)
library(FinTS)
library(rugarch)
library(dynlm) #for using lags in the model
library(vars) #for using vars
library(nlWaldTest) #for testing non linear wald test
library(lmtest) #for BP Test
library(broom) #for table presentations
library(vars) #for robust standard errors
library(sandwich)
library(knitr)
library(tsbox)
library(vrtest)
library(ggpubr)
library(stargazer)



#1. Read data
covid<- read_csv("COVIDcases.csv")
tweets<- read_csv("twitter.csv")
gdelt<- read_csv("gdelt.csv")

# Create date object

poscovid<-covid[(covid$month>= "2020-02-29" & covid$month < "2020-05-29"),]
posgdelt<-gdelt[(gdelt$Date>= "2020-02-29" & gdelt$Date < "2020-05-29"),]
postweet<-tweets[(tweets$month>= "2020-02-29" & tweets$month < "2020-05-29"),]



poscovid$month <- as.Date(poscovid$month, "%Y-%m-%d")
posgdelt$Date <- as.Date(posgdelt$Date, "%Y-%m-%d")
postweet$month <- as.Date(postweet$month, "%Y-%m-%d")

mycovid<-ts(poscovid$Canada, frequency = 7)
mygdelt<-ts(posgdelt$Canada, frequency = 7)
mytweet<-ts(postweet$Canada, frequency = 7)



par(mfrow=c(3,1))
plot(mycovid)
plot(mygdelt)
plot(mytweet)
plot(cbind(mycovid,mygdelt,mytweet))



decompco <- stl(log(mycovid), s.window="periodic")
plot(decompco)
ap.saco <- exp(seasadj(decompco))
autoplot(cbind(mycovid, SeasonallyAdjusted=ap.saco)) +
  xlab("Date") + ylab("Number of Covid cases")
mytscovid<-ap.saco

decompgd <- stl(log(mygdelt), s.window="periodic")
plot(decompgd)
ap.sagd <- exp(seasadj(decompgd))
autoplot(cbind(mygdelt, SeasonallyAdjusted=ap.sagd)) +
  xlab("Date") + ylab("Number of GDELT articles")
mytsgdelt<-ap.sagd

decomptw <- stl(log(mytweet), s.window="periodic")
plot(decomptw)
ap.satw <- exp(seasadj(decomptw))
autoplot(cbind(mytweet, SeasonallyAdjusted=ap.satw)) +
  xlab("Date") + ylab("Number of Tweets")
mytstweet<-ap.satw


# now to transform vector
#trans.vector = BoxCox(mycovid, lambda)

lambdacovid <- BoxCox.lambda(mytscovid)
lambdagdelt <- BoxCox.lambda(mytsgdelt)
lambdatweet <- BoxCox.lambda(mytstweet)

#pvalue of adtest should be less than 0.05
mycovid2 <- (mytscovid**(lambdacovid))
mycovid3 <- diff(mycovid2, lag=1)
adf.test(mycovid3)
ggtsdisplay(mycovid3)

mygdelt2 <- (mytsgdelt**(lambdagdelt))
mygdelt3 <- diff(mygdelt2, lag=1)
adf.test(mygdelt3)
ggtsdisplay(mygdelt3)


mytweet2 <- (mytstweet**(lambdatweet))
mytweet3 <- diff(mytweet2, lag=1)
adf.test(mytweet3)
#kpss.test(mytweet3)
ggtsdisplay(mytweet3)

plot(cbind(mycovid3,mygdelt3,mytweet3)) 


#choose lag given by AIC
v1 <- cbind(mycovid3, mygdelt3)

lagselect <- VARselect(v1, lag.max = 7, type = "const")
lagselect$selection #check p that meets normality, serial and stability
lagselect$criteria

#?VARselect


#Model1 <- VAR(v1, p = 8, type = "const", season = NULL, exog = NULL) 
Model1 <- VAR(v1, p = 6, type = "const", season = 7, exog = NULL) 
summary(Model1)
stargazer(Model1[["varresult"]], type='text')

#View(Model1)
#stability test two
roots(Model1, modulus = TRUE)

##The null hypothesis is no serial correlation,
#pvalue > than 0.05 to indicate lack of serial correlation
Serial1 <- serial.test(Model1, lags.pt = 40, type = "PT.asymptotic")
Serial1

plot(Serial1, names = "diffcovid")
plot(Serial1, names = "diffgdelt")


#The null hypothesis is no arch effects,
#test for heteroschedasticity , pvalue of > 0.05 indicate absence of heteroschedasticity

Arch1 <- arch.test(Model1, lags.multi = 42, multivariate.only = TRUE)
Arch1

#not normal, if pvalue<0.05, the model residuals are not normal
Norm1 <- normality.test(Model1, multivariate.only = TRUE)
Norm1

names(Model1)
#Plot of residuals

acf(residuals(Model1)[,1]) #covid
acf(residuals(Model1)[,2]) #gdelt

# test for the structural break in the residuals we can apply a CUSUM test.
#The stability test is seen if at any point in the graph, the sum goes out of 
#the red critical bounds, then a structural break at that point was seen.
Model1cusum <- stability(Model1, type = "OLS-CUSUM")
plot(Model1cusum)


resid <- residuals(Model1)
par(mfrow = c(1, 1))
plot.ts(resid[, 1])


xcroscovid_gdelt<-acf(residuals(Model1))
xcroscovid_gdelt$snames
acf(residuals(Model1)[,1])#acf for Covid, variable 1
acf(residuals(Model1)[,2])##acf for GDELT, variable 2 
acf(residuals(Model1)[,1:1])#acf plot 1 
acf(residuals(Model1)[,1:2])# acf plot 1,2, move clockwise
acf(residuals(Model1)[,2:1])# acf plot 3, move clockwise
acf(residuals(Model1)[,2:2])# plot 4, move clockwise

plot(xcroscovid_gdelt) # plots all acf residuals for var 1, 1&2, 2&1, 2

plot(xcroscovid_gdelt, xlim=c(0,20), ylim=c(-0.3,0.3),
     main="(a) CCF Plot Covid vs GDELT-Germany", ylab="Correlation",
     xlab="Time(days)")

xcroscovid_gdeltccf<-ccf(residuals(Model1)[,1],residuals(Model1)[,2],
                         xlim=c(-20,20), ylim=c(-0.3,0.3),
                         main="(a) CCF Plot Covid vs GDELT-Canada", ylab="Correlation",
                         xlab="Time(days)")

#twitter
v2 <- cbind(mycovid3, mytweet3)

lagselect <- VARselect(v2, lag.max = 7, type = "const")
lagselect$selection

#Model1 <- VAR(v1, p = 8, type = "const", season = NULL, exog = NULL) 
Model2 <- VAR(v2, p =5, type = "const", season = 7, exog = NULL) 
summary(Model2)

##The null hypothesis is no serial correlation,
#pvalue > than 0.05 to indicate lack of serial correlation
Serial2 <- serial.test(Model2, lags.pt = 60, type = "PT.asymptotic")
Serial2

#The null hypothesis is no arch effects,
#test for heteroschedasticity , pvalue of > 0.05 indicate absence of heteroschedasticity

Arch2 <- arch.test(Model2, lags.multi = 60, multivariate.only = TRUE)
Arch2



Model2cusum <- stability(Model2, type = "OLS-CUSUM")
plot(Model2cusum)
#Plot of residuals
acf(residuals(Model2)[,1]) #covid
acf(residuals(Model2)[,2]) #tweet

xcroscovid_tweet<-acf(residuals(Model2))
plot(xcroscovid_tweet)
plot(xcroscovid_tweet, xlim=c(0,20), ylim=c(-0.3,0.3),
     main="(b) CCF Plot Covid vs Twitter-Canada", ylab="Correlation",
     xlab="Time(days)")
#plot(Model2.fit)

xcroscovid_tweetccf<-ccf(residuals(Model2)[,1],residuals(Model2)[,2],
                         xlim=c(-20,20), ylim=c(-0.35,0.35),
                         main="(b) CCF Plot Covid vs Twitter-Canada", ylab="Correlation",
                         xlab="Time(days)")



tiff("Canadavar.tiff", units="in", width=9.7, height=4.0,compression="lzw", res=300)
par(mfrow=c(1,2), mar = c(4.6, 3.6, 1.6, 1.6))
#par(mfrow=c(1,2))
xcroscovid_gdeltccf<-ccf(residuals(Model1)[,1],residuals(Model1)[,2], 
                         xlim=c(-20,20), ylim=c(-0.35,0.35),
                         cex.axis=1.0, las=0.2,
                         main="",
                         ylab="",
                         xlab="",xaxs="i", yaxs="i" )
title (xlab="Time (weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(a) CCF Plot COVID-19 vs GDELT-Canada"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)


xcroscovid_tweetccf<-ccf(residuals(Model2)[,1],residuals(Model2)[,2],
                         xlim=c(-20,20), ylim=c(-0.35,0.35),
                         cex.axis=1.0, las=0.2,
                         main="", 
                         ylab="",
                         xlab="",)
title (xlab="Time(weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(b) CCF Plot COVID-19 vs Twitter-Canada"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)

dev.off() 

