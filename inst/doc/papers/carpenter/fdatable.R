#rm(list=ls())
load("matchfdaLin.RData")
sm <- summary(m.out)
tb <- cbind(sm$sum.matched[,3], sm$reduction[,1],
#            sm$sum.matched[,4], sm$reduction[,2],
            sm$sum.matched[,5], sm$reduction[,3],
            sm$sum.matched[,6], sm$reduction[,4])
colnames(tb) <- c(c(colnames(sm$sum.matched)[3],
                    colnames(sm$reduction)[1]),
#                  c(colnames(sm$sum.matched)[4],
#                    colnames(sm$reduction)[2]),
                  c(colnames(sm$sum.matched)[5],
                    colnames(sm$reduction)[3]),
                  c(colnames(sm$sum.matched)[6],
                    colnames(sm$reduction)[4]),
                  )
rownames(tb) <- rownames(sm$sum.all)
tb <- tb[c("distance", "prevgenx", "lethal", "deathrt1", "acutediz",
           "hosp01", "hospdisc", "femdiz01", "mandiz01", "peddiz01",
           "orphdum", "natreg", #"I(natreg^2)",
           "vandavg3", "wpnoavg3", "condavg3",
           "orderent", "stafcder"),]

rownames(tb)[rownames(tb)=="distance"] <- "Estimated propensity score"
rownames(tb)[rownames(tb)=="orderent"] <- "Order of disease market entry"
rownames(tb)[rownames(tb)=="prevgenx"] <- "Incidence of primary indication"
rownames(tb)[rownames(tb)=="lethal"] <- "Primary indication is lethal condition"
rownames(tb)[rownames(tb)=="deathrt1"] <- "Death rate, primary indication"
rownames(tb)[rownames(tb)=="hosp01"] <- "Primary ind. results in hospitalization"
rownames(tb)[rownames(tb)=="hospdisc"] <- "Hospitalizations associated with indication"
rownames(tb)[rownames(tb)=="hhosleng"] <- "Avg. length of hospitalization"
rownames(tb)[rownames(tb)=="femdiz01"] <- "Disease mainly affects women"
rownames(tb)[rownames(tb)=="mandiz01"] <- "Disease mainly affects men"
rownames(tb)[rownames(tb)=="peddiz01"] <- "Disease mainly affects children"
rownames(tb)[rownames(tb)=="acutediz"] <- "Primary indication is acute condition"
rownames(tb)[rownames(tb)=="natreg"] <- "National and regional groups"
#rownames(tb)[rownames(tb)=="I(natreg^2)"] <- "Square of national and regional groups"
rownames(tb)[rownames(tb)=="wpnoavg3"] <- "Washington post disease stories"
rownames(tb)[rownames(tb)=="vandavg3"] <- "Nightly TV news disease stories"
rownames(tb)[rownames(tb)=="condavg3"] <- "Days of Congressional hearings"
rownames(tb)[rownames(tb)=="stafcder"] <- "CDER staff"
rownames(tb)[rownames(tb)=="orphdum"] <- "Orphan Drug"

library(xtable)
print(xtable(tb))
