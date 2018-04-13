# Exploratory analyses of population-level effects

library(data.table)
library(ggplot2)


## Main effects analysis ##

dat = fread("../results/effect_size_data/model1maineffects_dailyfactor.csv")

unqpigs = dat[, list(pigID, study, cropuser)]
unqpigs = unqpigs[!duplicated(unqpigs)]
unqpigs[, list(users=sum(cropuser), total=length(cropuser)), by=study]


dat$beta[abs(dat$beta) > 4] = 0 # set large crop effects to 0

# Exploratory plots of crop users vs non-cropusers
meandat = dat[, list(betamean=mean(beta)), by=list(coef, study, cropuser)]
tplot = ggplot(meandat) + geom_boxplot(aes(x=coef, y=betamean, fill=cropuser)) + 
                      geom_hline(aes(yintercept=0)) + theme_bw() + 
                      xlab("Coefficient Name") + ylab("Effect size") + 
                      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
                      facet_wrap(~coef, scales="free")
ggsave("temp.pdf", width=20, height=20)



## Daily effects analysis ##

ddat = fread("../results/effect_size_data/model2dailyeffects_dailyfactor.csv")
ddat$hour = factor(ddat$hour, levels=c("morning", "midday", "evening"))


meanpig = ddat[, list(betamean = mean(beta)), by=list(coef, hour, study, cropuser)]
tplot = ggplot(meanpig) + geom_boxplot(aes(x=hour, y=betamean, fill=cropuser)) + 
                      geom_hline(aes(yintercept=0)) + theme_bw() + 
                      xlab("Coefficient Name") + ylab("Effect size") + 
                      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
                      facet_wrap(~coef, scales='free') #+ guides(color=FALSE)
ggsave("temp2.pdf", width=12, height=12)


## Seasonal effects analysis ##




