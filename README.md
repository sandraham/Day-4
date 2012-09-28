# Assignment 4


### Load libraries and data

library("TraMineR")

data(biofam)

cohort <- cut(biofam$birthyr, breaks=c(1900,1930,1940,1950,1960), 
              labels=c("1900-1929","1930-1939","1940-1949","1950-1959"), right=FALSE)

biofam2 <- cbind(biofam, cohort)


## 1. Build a table of characteristics to use for distributional plots
tab <- data.frame(seqlength(biofam2.seq), seqtransn(biofam2.seq), seqsubsn(biofam2.seq),
                   seqient(biofam2.seq), seqST(biofam2.seq), seqici(biofam2.seq))


# 2. Get distribution statistics for the computed longitudinal characteristics

summary(tab, digits=3)


# 3. Histograms of longitudinal characteristics

par(mfrow=c(3,2))
hist(tab$Trans.,  col="LightGreen", main="# Transitions")
hist(tab$Subseq.,  col="DarkGreen", main="# Subsequences")
hist(tab$Entropy,  col="Lavender", main="Entropy")
hist(tab$Turbulence,  col="LightBlue", main="Turbulence")
hist(tab$C,  col="Blue", main="Complexity Index")



# 4. List distinct successive states (DSS)
seqdss(biofam2.seq)[1995:2000, ]

# Duration of successive states
seqdur(biofam2.seq)[1995:2000, ]


# 5. Mean and variance of time spent in successive states
dur <- seqdur(biofam2.seq)

dur.mean <- apply(dur, 1, mean,  na.rm=TRUE)
dur.var <- apply(dur, 1, var,  na.rm=TRUE)

summary(dur.mean)
summary(dur.var)


# 6. Scatterplot matrix of entropy, turbulence, and the complexity index

plot(tab[, 2:4], main="Scatterplot Matrix of Longitudinal Characteristics")


# 7. Boxplots by birth cohort

boxplot(tab$C ~ biofam2$cohort, col="LightBlue", main="Complexity Index")


# 8. Regressions for complexity on covariates

lm.C <- lm(tab$C ~ cohort + sex + plingu02, data=biofam2)

summary(lm.C)

The mean complexity index for French-speaking men born between 1900 and 1929 is 0.15,
indicating fairly simple family sequences with fewer transitions and lower entropy.
  Later cohorts had significantly higher complexity indices, by 0.04 for those born
in the 1930's, 0.05 for those born in the 1940's, and 0.06 for those born in the 
1950's (all p=0.001).
   
Women had complexity indices that average 0.01 more than men (p=0.01) and 
Italian-speakers had indices 0.05 lower than those who interviewed in French 
or German (p=0.001).

