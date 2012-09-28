> # Assignment 4
> 
> 
> ### Load libraries and data
> 
> library("TraMineR")
> 
> data(biofam)
> 
> cohort <- cut(biofam$birthyr, breaks=c(1900,1930,1940,1950,1960), 
+               labels=c("1900-1929","1930-1939","1940-1949","1950-1959"), right=FALSE)
> 
> biofam2 <- cbind(biofam, cohort)
> 
> 
> ## 1. Build a table of characteristics to use for distributional plots
> tab <- data.frame(seqlength(biofam2.seq), seqtransn(biofam2.seq), seqsubsn(biofam2.seq),
+                    seqient(biofam2.seq), seqST(biofam2.seq), seqici(biofam2.seq))
 [>] computing entropy for 2000 sequences ...
 [>] computing state distribution for 2000 sequences ...
 [>] extracting symbols and durations ...
 [>] computing turbulence for 2000 sequence(s) ...
 [>] computing complexity index for 2000 sequences ...
> 
> 
> # 2. Get distribution statistics for the computed longitudinal characteristics
> 
> summary(tab, digits=3)
     Length       Trans.        Subseq.         Entropy        Turbulence         C        
 Min.   :16   Min.   :0.00   Min.   : 2.00   Min.   :0.000   Min.   :1.00   Min.   :0.000  
 1st Qu.:16   1st Qu.:1.00   1st Qu.: 4.00   1st Qu.:0.299   1st Qu.:3.69   1st Qu.:0.141  
 Median :16   Median :1.00   Median : 4.00   Median :0.333   Median :5.06   Median :0.149  
 Mean   :16   Mean   :1.56   Mean   : 7.05   Mean   :0.355   Mean   :4.80   Mean   :0.191  
 3rd Qu.:16   3rd Qu.:2.00   3rd Qu.: 8.00   3rd Qu.:0.473   3rd Qu.:6.22   3rd Qu.:0.251  
 Max.   :16   Max.   :4.00   Max.   :32.00   Max.   :0.703   Max.   :8.81   Max.   :0.433  
> 
> 
> # 3. Histograms of longitudinal characteristics
> 
> par(mfrow=c(3,2))
> hist(tab$Trans.,  col="LightGreen", main="# Transitions")
> hist(tab$Subseq.,  col="DarkGreen", main="# Subsequences")
> hist(tab$Entropy,  col="Lavender", main="Entropy")
> hist(tab$Turbulence,  col="LightBlue", main="Turbulence")
> hist(tab$C,  col="Blue", main="Complexity Index")
> 
> 
> 
> # 4. List distinct successive states (DSS)
> seqdss(biofam2.seq)[1995:2000, ]
     Sequence  
59   P-M       
629  P-L-LMC   
2297 P-L-LM-LMC
775  P         
2522 P-M       
719  P-LMC     
> 
> # Duration of successive states
> seqdur(biofam2.seq)[1995:2000, ]
     DUR1 DUR2 DUR3 DUR4 DUR5
59     13    3   NA   NA   NA
629     6    3    7   NA   NA
2297    2    6    4    4   NA
775    16   NA   NA   NA   NA
2522    3   13   NA   NA   NA
719     6   10   NA   NA   NA
> 
> 
> # 5. Mean and variance of time spent in successive states
> dur <- seqdur(biofam2.seq)
> 
> dur.mean <- apply(dur, 1, mean,  na.rm=TRUE)
> dur.var <- apply(dur, 1, var,  na.rm=TRUE)
> 
> summary(dur.mean)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  3.200   5.333   8.000   7.129   8.000  16.000 
> summary(dur.var)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   4.667  11.330  18.820  20.330  98.000     162 
> 
> 
> # 6. Scatterplot matrix of entropy, turbulence, and the complexity index
> 
> plot(tab[, 2:4], main="Scatterplot Matrix of Longitudinal Characteristics")
> 
> 
> # 7. Boxplots by birth cohort
> 
> boxplot(tab$C ~ biofam2$cohort, col="LightBlue", main="Complexity Index")
> 
> 
> # 8. Regressions for complexity on covariates
> 
> lm.C <- lm(tab$C ~ cohort + sex + plingu02, data=biofam2)
> summary(lm.C)

Call:
lm(formula = tab$C ~ cohort + sex + plingu02, data = biofam2)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.218045 -0.061923  0.002679  0.054963  0.230796 

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)      0.146918   0.007914  18.563  < 2e-16 ***
cohort1930-1939  0.039250   0.007964   4.929 9.12e-07 ***
cohort1940-1949  0.052539   0.007557   6.952 5.16e-12 ***
cohort1950-1959  0.057482   0.007554   7.609 4.62e-14 ***
sexwoman         0.013645   0.004532   3.011  0.00265 ** 
plingu02german  -0.003482   0.005134  -0.678  0.49765    
plingu02italian -0.047038   0.010980  -4.284 1.94e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

Residual standard error: 0.09127 on 1640 degrees of freedom
  (353 observations deleted due to missingness)
Multiple R-squared: 0.05307,  Adjusted R-squared: 0.04961 
F-statistic: 15.32 on 6 and 1640 DF,  p-value: < 2.2e-16 

> 
>   ##### The mean complexity index for French-speaking men born between 1900 and 1929 is 0.15,
>   #       indicating fairly simple family sequences with fewer transitions and lower entropy.
>   #
>   #       Later cohorts had significantly higher complexity indices, by 0.04 for those born
>   #       in the 1930's, 0.05 for those born in the 1940's, and 0.06 for those born in the 
>   #       1950's (all p=0.001).
>   # 
>   #       Women had complexity indices that average 0.01 more than men (p=0.01) and 
>   #       Italian-speakers had indices 0.05 lower than those who interviewed in French 
>   #       or German (p=0.001).
>   ############################