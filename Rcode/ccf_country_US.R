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

plot(poscovid$month,poscovid$`United States`, xaxt = "n", type = "l")
axis(1, poscovid$month, format(poscovid$month, "%Y-%m-%d"), cex.axis = .7)

# First we have to change the date to POSIXct
poscovid$month <- strptime(poscovid$month, "%Y-%m-%d" )
poscovid$month<- as.POSIXct(poscovid$month)
#poscovid$`United States` <- as.numeric(poscovid$`United States`)

mycovid<-ts(poscovid$`United States`, frequency = 7)
mygdelt<-ts(posgdelt$`United States`, frequency = 7)
mytweet<-ts(postweet$`United States`, frequency = 7)





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


mygdelt2 <- (mytsgdelt**(1/4))
mygdelt3 <- mygdelt2
adf.test(mygdelt3)
ggtsdisplay(mygdelt3)


mytweet2 <- (mytstweet**(1/2))
mytweet3 <- mytweet2
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


tiff("Fig4.tiff", units="in", width=6.3, height=6.0,compression="lzw", res=300)

par(mfrow=c(2,1),  mai = c(0.8, 1.3, 0.2, 0.5))
plot(mycovid3, xaxt='n', xlab=" ", 
     ylab="Stationary data \nCOVID-19 \ncases") 
text="(a)"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 0.8, col = NA, font = 1)
plot(poscovid$month, mycovid, xlab="Date",
     ylab="Original data \nCOVID-19 \ncases", type='l') 
text="(b) "
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 0.8, col = NA, font = 1)

#dev.off()  
graphics.off() 



fitdiffcovid<-auto.arima(mycovid3, stepwise=FALSE, approximation=FALSE,
                         seasonal = TRUE)
tsdiag(fitdiffcovid)

tiff("Fig5.tiff", units="in", width=6.5, height=4.75, compression="lzw", res=300)
acf(fitdiffcovid$residuals, main="ARIMA model for COVID-19 residuals", ylab="Correlation",
    xlab="Time (weeks)", cex.axis=0.8)
graphics.off()

checkresiduals(fitdiffcovid, main = "",
               ylab="a",
               xlab="Time(days)") #pvalue should be > 0.05 

par(mfrow=c(3,1))
plot(fitdiffcovid$residuals)
acf(fitdiffcovid$residuals)
pacf(fitdiffcovid$residuals)

#Step 3- store residuals to a variable
xresdiffcovid<-fitdiffcovid$residuals


#gdelt di
ydiffgdeltfitted<-fitted(Arima(mygdelt3,model=fitdiffcovid))
yresdiffgdelt <- mygdelt3- fitted(Arima(mygdelt3,model=fitdiffcovid))


xcroscovid_gdelt<-ccf(xresdiffcovid,yresdiffgdelt,na.action=na.omit, 
                      xlim=c(-3,3), ylim=c(-0.3,0.3), main="CCF Plot Covid vs GDELT-U.S", ylab="Correlation",
                      xlab="Time(weeks)")
text(1.5, 0.1, "1", col='red', cex=2)
text(1.5, -0.1, "2", col='red', cex=2)
text(-1.5, -0.1, "3", col='red', cex=2)
text(-1.5, 0.1, "4", col='red', cex=2)
abline(v=0, col="blue", lwd=3, lty=1)
abline(h=0, col="blue", lwd=3, lty=1)
#?ccf

xcroscovid_gdeltdf<-data.frame(xcroscovid_gdelt$lag,xcroscovid_gdelt$acf)
xcroscovid_gdeltresdf<-data.frame(lag=xcroscovid_gdelt$lag*7,acfvalues=xcroscovid_gdelt$acf)

orderxcroscovid_gdelt = xcroscovid_gdeltresdf[order(-xcroscovid_gdeltresdf$acf),]
top2xcrosscovid_gdelt = head(orderxcroscovid_gdelt, 2)
top2xcrosscovid_gdelt 



#Tweet

ydifftweetfitted<-fitted(Arima(mytweet3,model=fitdiffcovid))
yresdifftweet <- mytweet3- fitted(Arima(mytweet3,model=fitdiffcovid))

xcroscovid_tweet<-ccf(xresdiffcovid,yresdifftweet,na.action=na.omit,
                      xlim=c(-3,3), ylim=c(-0.3,0.3), main="CCF Plot Covid vs Twitter-U.S", 
                      ylab="Correlation",
                      xlab="Time(weeks")
text(1.5, 0.1, "1", col='red', cex=2)
text(1.5, -0.1, "2", col='red', cex=2)
text(-1.5, -0.1, "3", col='red', cex=2)
text(-1.5, 0.1, "4", col='red', cex=2)
abline(v=0, col="blue", lwd=3, lty=2)
abline(h=0, col="blue", lwd=3, lty=2)

xcrosresultcovid_tweet<-data.frame(xcroscovid_tweet$lag,xcroscovid_tweet$acf)
xcrosresultcovid_tweetdf<-data.frame(lag=xcroscovid_tweet$lag*7,acfvalues=xcroscovid_tweet$acf)

orderxcroscovid_tweet = xcrosresultcovid_tweetdf[order(-xcroscovid_tweet$acf),]
top2xcrosscovid_tweet = head(orderxcroscovid_tweet, 2)
top2xcrosscovid_tweet 

#?par
#plots

tiff("US.tiff", units="in", width=9.5, height=4.0,compression="lzw", res=300)
par(mfrow=c(1,2), mar = c(4.6, 3.6, 1.6, 1.6))

xcroscovid_gdelt<-ccf(xresdiffcovid,yresdiffgdelt,na.action=na.omit, 
                      xlim=c(-3,3), ylim=c(-0.35,0.35),
                      cex.axis=1.0, las=0.2,
                      main="",
                      ylab="",
                      xlab="",xaxs="i", yaxs="i" )
title (xlab="Time (weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(a) CCF Plot COVID-19 vs GDELT-U.S"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)
text(1.5, 0.1, "1", col='red', cex=1.5)
text(1.5, -0.1, "2", col='red', cex=1.5)
text(-1.5, -0.1, "3", col='red', cex=1.5)
text(-1.5, 0.1, "4", col='red', cex=1.5)
abline(v=0, col="blue", lwd=3, lty=1)
abline(h=0, col="blue", lwd=3, lty=1)


xcroscovid_tweet<-ccf(xresdiffcovid,yresdifftweet,na.action=na.omit,
                      xlim=c(-3,3), ylim=c(-0.35,0.35),
                      cex.axis=1.0, las=0.2,
                      main="", 
                      ylab="",
                      xlab="",)
title (xlab="Time(weeks)", cex.lab=1.2, line=2.0)
title(ylab="Correlation", cex.lab=1.2, line=2.0)
text="(b) CCF Plot COVID-19 vs Twitter-U.S"
mtext(text, side = 3, line = 0, outer = FALSE, at =NA, adj =0, padj = 0, cex = 1.2, col = NA, font = 1)

dev.off() 