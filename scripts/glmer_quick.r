library(ape)

alllatlon <- read.table(file=file.path(getwd(),'outputs','quick_alr','latlon_north_america_phy.txt'),sep='\t', row.names=NULL)

phy_time <- read.tree(file.path(getwd(),'outputs','quick_alr','phy_filtered.newick'))
phy_divergence <- read.tree(file.path(getwd(),'datasets/raw/NextStrain/nextstrain_ncov_north-america_tree.nwk'))

phy_divergence <- drop.tip(phy_divergence, phy_divergence$tip.label[!phy_divergence$tip.label %in% phy_time$tip.label])

rates <- phy_divergence$edge.length / phy_time$edge.length

alllatlon$rate <- 0
alllatlon$rate[1:length(phy_time$tip.label)] <- rates[1:length(phy_time$tip.label)]
alllatlon$rate[(length(phy_time$tip.label)+2):nrow(alllatlon)] <- rates[(length(phy_time$tip.label)+1):(nrow(alllatlon)-1)] 

write.table(alllatlon,file=file.path(getwd(),'outputs','quick_alr','latlon_north_america_phy_divergence_rate.txt'),sep='\t',quote=FALSE,row.names=FALSE)

pdf(file.path(getwd(),'outputs','quick_alr','log_rates.pdf'))
hist(log(rates),100,main='Log rate of evolution')
dev.off()

pdf(file.path(getwd(),'outputs','quick_alr','rates.pdf'))
hist(rates,100,main='Rate of evolution')
dev.off()

#preds <- read.table(file.path(getwd(),'outputs/gp/gp_predictions_NOAA.txt'),sep='\t', row.names=NULL, header=TRUE)
#cov_title <- 'NOAA Temperature Anomalies Jan 2017'
preds <- read.table(file.path(getwd(),'outputs/gp/gp_predictions_NOAA.txt'),sep='\t', row.names=NULL, header=TRUE)
cov_title <- 'CDC Particulate Matter Concentration'
title_short <- 'pollution'
covariate <- apply(alllatlon[-(length(phy_divergence$edge.length)+1),c(3,2)],1, function(x) preds[apply(preds[,2:3],1,function(y) all(y == x)), 1])

#covariate <- rnorm(length(phy_time$edge.length))

dat <- as.data.frame(cbind(cbind(phy_divergence$edge.length, log(phy_time$edge.length)), covariate))
colnames(dat) <- c('div','time','co')

res <- glmer(div ~ co + (1|time), data = dat, family='poisson')
summary(res)

pdf(file.path(getwd(),'outputs','quick_alr',paste0('mutations_vs_', title_short, '.pdf')))
plot(dat$co,dat$div, xlab = cov_title, ylab = 'Number of mutations')
dev.off()

pdf(file.path(getwd(),'outputs','quick_alr',paste0('log_mutations_vs_', title_short, '.pdf')))
plot(dat$co,log(dat$div), xlab = cov_title, ylab = 'Log of number of mutations')
dev.off()
                                                                                                             
pdf(file.path(getwd(),'outputs','quick_alr','mutations_vs_time.pdf'))
plot(exp(dat$time),dat$div, xlab = 'Time', ylab = 'Number of mutations')
dev.off()

