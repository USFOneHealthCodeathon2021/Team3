library(ape)
library(cmdstanr)
library(rstan)

wd <- getwd()
output_prefix <- file.path(wd,'outputs/alr')
model_dir <- file.path(wd, 'models')
model_name <- 'alr'

sampling_command <- paste(paste0('./', model_name),
                paste0('data file=', file.path(output_prefix, 'data.json')),
                #paste0('init=', file.path(output_prefix, 'inits.json')),
                'output',
                paste0('file=', file.path(output_prefix, 'samples.txt')),
                paste0('refresh=', 100),
                'method=sample',
                'adapt delta=0.95',
                'algorithm=hmc',
                #'stepsize=0.01',
                'engine=nuts',
                'max_depth=10',
                'num_warmup=1000',
                'num_samples=200',
                sep=' ')
      
phy <- read.tree(file.path(getwd(),'outputs','quick_alr','phy_filtered.newick'))
latlon_unique <- read.table(file.path(getwd(),'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',header=TRUE)
dat <- read.table(file.path(getwd(),'datasets','curated','NextStrain','nextstrain_ncov_north-america_metadata.tsv'), header=T, sep='\t', quote='')
dat <- dat[dat$Location != '',]

latlon_tips <- t(sapply(phy$tip.label,function(x) latlon_unique[latlon_unique$query==dat[dat$Strain==x,'Location'],c('lon','lat')]))


phy <- drop.tip(phy,phy$tip.label[!phy$tip.label %in% rownames(latlon_tips)])
latlon <- matrix(as.numeric(latlon_tips[phy$tip.label,]),ncol=2)
latlonrad <- DescTools:::DegToRad(latlon)
      
N <- nrow(latlon)
NI <- phy$Nnode

tipdist <- geosphere:::distm(latlon, fun = function(x,y) geosphere:::distGeo(x,y,a=1))
ratevec <- as.vector((tipdist/cophenetic(phy))[lower.tri(tipdist)])
                             
sigma_prior <- mean(ratevec[!is.infinite(ratevec) & !is.na(ratevec)]) * 10
                             
                             
data <- list(N = N,
             NI = NI,
             phi_tips = latlonrad[,1],
             theta_tips = latlonrad[,2],
             time = phy$edge.length + 1e-9,
             self = phy$edge[,2],
             ancestor = phy$edge[,1],
             sigma_prior = sigma_prior)

#init <- list()

#write_stan_json(init, file.path(output_prefix, 'inits.json'))
write_stan_json(data, file.path(output_prefix, 'data.json'))

setwd(cmdstan_path())
system(paste0('make ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_command)
print(date())
system(sampling_command)

setwd(wd)
stanfit <- read_stan_csv(file.path(output_prefix, 'samples.txt'))
check_hmc_diagnostics(stanfit)
                             
locations_cartesian <- apply(extract(stanfit,pars='loc_anc')[[1]],c(2,3),mean)
locations_radians <- t(apply(locations_cartesian, 1, function(x) c(phi=atan2(x[2],x[1]),theta=atan2(sqrt(sum(x[1:2]^2)),x[3]))))
locations_degrees <- locations_radians / pi * 180
locations_degrees[,2] <- locations_degrees[,2] - 90
colnames(locations_degrees) <- c('longitude','latitude')
                             
sigma <- mean(extract(stanfit,pars='sigma_raw')[[1]])
