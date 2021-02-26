#install.packages('geosphere')
#install.packages('DescTools')
library(cmdstanr)
library(rstan)

wd <- getwd()
output_prefix <- file.path(wd,'outputs/gp')
model_dir <- file.path(wd, 'models')
model_name <- 'gp_continuous'
dataset_name <- 'pollution'

sampling_command <- paste(paste0('./', model_name),
                paste0('data file=', file.path(output_prefix, 'data.json')),
                paste0('init=', file.path(output_prefix, 'inits.json')),
                #'init=0',
                'output',
                paste0('file=', file.path(output_prefix, paste0('samples_',dataset_name,'.txt'))),
                paste0('refresh=', 100),
                'method=sample algorithm=hmc',
                #'stepsize=0.01',
                'engine=nuts',
                'max_depth=7',
                'num_warmup=500',
                'num_samples=200',
                sep=' ')
                
                
testdat <- read.table(file.path(wd, 'datasets/curated/NOAAGlobalTemp/testdat.txt'),sep='\t',header=T)
pollutiondat <- read.csv(file.path(wd, 'datasets/curated/CDC_PM2.5_Concentrations/ave-county-level-pm2.5-2014.csv'),header=T)
library('tmaptools')
pollutiondat$combined <- apply(pollutiondat,1,function(x) paste(x[2],x[1],sep=', '))
#latlon <- geocode_OSM(pollutiondat$combined)
#write.table(latlon,file=file.path(wd,'datasets','curated','CDC_PM2.5_Concentrations','latlon_ave-county-level-pm2.5-2014.txt'),sep='\t',quote=FALSE,row.names=FALSE)
latlon <- read.table(file.path(wd,'datasets','curated','CDC_PM2.5_Concentrations','latlon_ave-county-level-pm2.5-2014.txt'),sep='\t',header=TRUE)
pollutiondat <- pollutiondat[pollutiondat$combined %in% latlon$query,]
N <- nrow(latlon)

alllatlon <- read.table(file.path(getwd(),'outputs','quick_alr','latlon_north_america_phy.txt'), sep='\t',quote='',header=TRUE,row.names=NULL)
predictLatLon <- unique(alllatlon[,c(3,2)])
N2 <- nrow(predictLatLon)

#N2 <- 100
#randomCartesian <- t(apply(matrix(rnorm(N2*3), nrow=N2), 1, function(x) x / sqrt(sum(x^2))))                   
#randomPolar <- t(apply(randomCartesian, 1, function(x) unlist(DescTools:::CartToSph(x[1],x[2],x[3]))))                       
#predictLatLon <- DescTools:::RadToDeg(randomPolar[,2:3])

colnames(predictLatLon) <- c('lat','lon')

d2 <- geosphere:::distm(rbind(latlon[,c(3,2)],predictLatLon[,c(2,1)]), fun = function(x,y) geosphere:::distGeo(x,y,a=1))

y_norm <- pollutiondat[,3] / sd(pollutiondat[,3])                       
                        
data <- list(N = N,
             d = d2[1:N,1:N],
             prior_scale = mean(d2[lower.tri(d2)]),
             N2 = N2,
             d2 = d2[(N+1):(N+N2),],
             y = y_norm)

init <- list(rho_raw=0.1,alpha=1,sigma=1)

write_stan_json(init, file.path(output_prefix, 'inits.json'))
write_stan_json(data, file.path(output_prefix, 'data.json'))

setwd(cmdstan_path())
system(paste0('make ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_command)
print(date())
system(sampling_command)

setwd(wd)
stanfit <- read_stan_csv(file.path(output_prefix, paste0('samples_',dataset_name,'.txt')))
check_hmc_diagnostics(stanfit)
summary(stanfit)
                        
y2 <- apply(extract(stanfit, pars='y2')[[1]],2,mean)
pollution_predictions <- cbind(y2,predictLatLon)
                       
write.table(pollution_predictions, file=file.path(output_prefix,'gp_predictions_pollution.txt'), sep='\t', quote=FALSE, row.names=FALSE)

                        
                        
                        
