library(TSA)
library(stats)
library(readr)
library(tidyverse)
library(tseries)
library(forecast)
library(fpp2)
library(gam)
library(lmtest)
library(xts)
#library(quantmod)
#1.Read data
covid<- read_csv("COVIDcases.csv")
tweets<- read_csv("twitter.csv")
gdelt<- read_csv("gdelt.csv")

#2. Create date object
#Date filter
poscovid<-covid[(covid$month>= "2020-02-29" & covid$month < "2020-05-29"),]
posgdelt<-gdelt[(gdelt$Date>= "2020-02-29" & gdelt$Date < "2020-05-29"),]
postweet<-tweets[(tweets$month>= "2020-02-29" & tweets$month < "2020-05-29"),]


poscovid$month <- as.Date(poscovid$month, "%Y-%m-%d")
posgdelt$Date <- as.Date(posgdelt$Date, "%Y-%m-%d")
postweet$month <- as.Date(postweet$month, "%Y-%m-%d")


# First we have to change the date to POSIXct
poscovid$month <- strptime(poscovid$month, "%Y-%m-%d" )
poscovid$month<- as.POSIXct(poscovid$month)
poscovid$Italy <- as.numeric(poscovid$Italy)

mycovid<-ts(poscovid$Italy, frequency = 7)
mygdelt<-ts(posgdelt$Italy, frequency = 7)
mytweet<-ts(postweet$Italy, frequency = 7)


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
mycovid2 <- (mytscovid**(1/256))
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
ggtsdisplay(mytweet3)

plot(cbind(mycovid3,mygdelt3,mytweet3))  

par(mfrow=c(3,2))
plot((mycovid)) 
plot(mycovid3) 
plot((mygdelt)) 
plot(mygdelt3) 
plot((mytweet)) 
plot(mytweet3) 


fitdiffcovid<-auto.arima(mycovid3, stepwise=FALSE, approximation=FALSE,
                         seasonal = TRUE)
tsdiag(fitdiffcovid)
checkresiduals(fitdiffcovid) #pvalue should be > 0.05

par(mfrow=c(3,1))
plot(fitdiffcovid$residuals)
acf(fitdiffcovid$residuals)
pacf(fitdiffcovid$residuals)

#Step 3- store residuals to a variable
xresdiffcovid<-fitdiffcovid$residuals

#gdelt
yresdiffgdelt <- mygdelt3- fitted(Arima(mygdelt3,model=fitdiffcovid))

xcroscovid_gdelt<-ccf(xresdiffcovid,yresdiffgdelt,na.action=na.omit, 
                      xlim=c(-3,3), ylim=c(-0.3,0.3), main="CCF Plot Covid vs GDELT-Italy", ylab="Correlation",
                      xlab="Time(weeks)")
#?ccf

xcroscovid_gdeltdf<-data.frame(xcroscovid_gdelt$lag,xcroscovid_gdelt$acf)
xcroscovid_gdeltresdf<-data.frame(lag=xcroscovid_gdelt$lag*7,acfvalues=xcroscovid_gdelt$acf)

orderxcroscovid_gdelt = xcroscovid_gdeltresdf[order(-xcroscovid_gdeltresdf$acf),]
top2xcrosscovid_gdelt = head(orderxcroscovid_gdelt, 2)
top2xcrosscovid_gdelt 



#Tweet


yresdifftweet <- mytweet3- fitted(Arima(mytweet3,model=fitdiffcovid))
xcroscovid_tweet<-ccf(xresdiffcovid,yresdifftweet,na.action=na.omit,
                      xlim=c(-3,3), ylim=c(-0.3,0.3), main="CCF Plot Covid vs Twitter-Italy", 
                      ylab="Correlation",
                      xlab="Time(weeks)")
xcrosresultcovid_tweet<-data.frame(xcroscovid_tweet$lag,xcroscovid_tweet$acf)
xcrosresultcovid_tweetdf<-data.frame(lag=xcroscovid_tweet$lag*7,acfvalues=xcroscovid_tweet$acf)

orderxcroscovid_tweet = xcrosresultcovid_tweetdf[order(-xcroscovid_tweet$acf),]
top2xcrosscovid_tweet = head(orderxcroscovid_tweet, 2)
top2xcrosscovid_tweet 

#?par
tiff("Italy_coded.tiff", units="in", width=9.7, height=4.0, res=300)
par(mfrow=c(1,2), mar = c(4.6, 3.6, 1.6, 1.6))
#par(mfrow=c(1,2))
xcroscovid_gdelt<-ccf(xresdiffcovid,yresdiffgdelt,na.action=na.omit, 
                      xlim=c(-3,3), ylim=c(-0.3,0.3),
                      cex.axis=1.0, las=0.2,
                      main="",
                      ylab="",
                      xlab="",xaxs="i", yaxs="i" )
title (xlab="Time (weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(g) CCF Plot COVID-19 vs GDELT-Italy"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)


xcroscovid_tweet<-ccf(xresdiffcovid,yresdifftweet,na.action=na.omit,
                      xlim=c(-3,3), ylim=c(-0.3,0.3),
                      cex.axis=1.0, las=0.2,
                      main="", 
                      ylab="",
                      xlab="",)
title (xlab="Time(weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(h) CCF Plot COVID-19 vs Twitter-Italy"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)

dev.off() 