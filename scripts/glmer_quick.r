library(ape)

alllatlon <- read.table(file=file.path(getwd(),'outputs','quick_alr','latlon_north_america_phy.txt'),sep='\t', row.names=NULL)

phy_time <- read.tree(file.path(getwd(),'outputs','quick_alr','phy_filtered.newick'))
phy_divergence <- read.tree(file.path(getwd(),'datasets/raw/NextStrain/nextstrain_ncov_north-america_tree.nwk'))

phy_divergence <- drop.tip(phy_divergence, phy_divergence$tip.label[!phy_divergence$tip.label %in% phy_time$tip.label])

rates <- phy_divergence$edge.length / phy_time$edge.length

pdf(file.path(getwd(),'outputs','quick_alr','log_rates.pdf'))
hist(log(rates),100,main='Log rate of evolution')
dev.off()

pdf(file.path(getwd(),'outputs','quick_alr','rates.pdf'))
hist(rates,100,main='Rate of evolution')
dev.off()

covariate <- rnorm(length(phy_time$edge.length))

dat <- as.data.frame(cbind(cbind(phy_divergence$edge.length, log(phy_time$edge.length)), covariate))
colnames(dat) <- c('div','time','co')

res <- glmer(div ~ co + (1|time), data = dat, family='poisson')
summary(res)

pdf(file.path(getwd(),'outputs','quick_alr','mutations_vs_covariate.pdf'))
plot(dat$co,dat$div, xlab = 'Covariate', ylab = 'Number of mutations')
dev.off()

pdf(file.path(getwd(),'outputs','quick_alr','mutations_vs_time.pdf'))
plot(exp(dat$time),dat$div, xlab = 'Time', ylab = 'Number of mutations')
dev.off()

