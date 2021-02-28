library(ape)
library(cmdstanr)
library(rstan)

wd <- getwd()
output_prefix <- file.path(wd,'outputs/alr')
model_dir <- file.path(wd, 'models')
model_name <- 'alr_noise_donut'

sampling_command <- paste(paste0('./', model_name),
                paste0('data file=', file.path(output_prefix, paste0(model_name, '_data.json'))),
                paste0('init=', file.path(output_prefix, paste0(model_name, '_inits.json'))),
                'output',
                paste0('file=', file.path(output_prefix, paste0(model_name, '_samples.txt'))),
                paste0('refresh=', 100),
                'method=sample',
                'adapt delta=0.8',
                'algorithm=hmc',
                #'stepsize=0.01',
                'engine=nuts',
                'max_depth=10',
                'num_warmup=1000',
                'num_samples=500',
                sep=' ')


phy <- collapse.singles(read.tree(file.path(wd,'datasets/raw/NextStrain/nextstrain_ncov_north-america_timetree.nwk')))
dat <- read.table(file.path(wd,'datasets/raw/NextStrain/nextstrain_ncov_north-america_metadata.tsv'),header=T, sep='\t', quote='')

dat$location_combined <- apply(dat[,c('Location','Admin.Division','Country')], 1, function(x) paste(x,collapse=', '))
#latlon <- tmaptools:::geocode_OSM(unique(dat$location_combined))
#latlon <- latlon[!(is.na(latlon$lat) | is.na(latlon$lon)),]
#write.table(latlon,file=file.path(wd,'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',quote=FALSE,row.names=FALSE)
latlon <- read.table(file.path(getwd(),'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',header=TRUE)

latlon$kappa <- apply(latlon[,c('lon_min','lon_max','lat_min','lat_max')], 1, function(x) {
  r <- geosphere:::distGeo(c(x[[1]],x[[3]]),c(x[[2]],x[[4]]), a=1) / 2
  return(qnorm(0.99,0,1/sqrt(r)))
})

latlon_tips <- t(sapply(phy$tip.label, function(x) latlon[latlon$query==dat[dat$Strain==x,'location_combined'],c('lon','lat','kappa')]))                       
#phy <- drop.tip(phy,phy$tip.label[!phy$tip.label %in% rownames(latlon_tips)])

#write.tree(phy, file=file.path(output_prefix, 'nextstrain_ncov_north-america_timetree_filteredGPS.newick'))  
phy <- read.tree(file.path(output_prefix, 'nextstrain_ncov_north-america_timetree_filteredGPS.newick'))


latlon_mat <- matrix(as.numeric(latlon_tips[phy$tip.label,]),ncol=3)
latlonrad <- latlon_mat[,1:2] / 180 * pi
latlonrad[,2] <- (pi/2) - latlonrad[,2] ## geographic degrees go from -90 to +90 starting at equator, math starts at top and goes 0-pi
                        
latlon_internal_guess <- t(sapply(unique(phy$edge[,1]), function(x) icosa:::surfacecentroid(matrix(as.numeric(latlon_tips[phangorn:::Descendants(phy,x,'tips')[[1]],1:2]),ncol=2))))   
latlon_internal_guess_rad <- DescTools:::DegToRad(latlon_internal_guess)
latlon_internal_guess_rad[,2] <- (pi/2) - latlon_internal_guess_rad[,2]
                                  
latlon_inits_rad <- rbind(latlonrad, latlon_internal_guess_rad)
latlon_inits_cartesian <- t(apply(latlon_inits_rad, 1, function(x) c(sin(x[2])*cos(x[1]), sin(x[2])*sin(x[1]), cos(x[2])))) 
latlon_inits_cartesian_noise <- latlon_inits_cartesian + matrix(rnorm(nrow(latlon_inits_cartesian)*3,0,1e-4),ncol=3)
latlon_inits_cartesian_noise <- t(apply(latlon_inits_cartesian_noise, 1, function(x) x / sqrt(sum(x^2))))                                   
      
N <- nrow(latlon_mat)
NI <- phy$Nnode

tipdist <- geosphere:::distm(latlon_mat[,1:2], fun = function(x,y) geosphere:::distGeo(x,y,a=1))

time_min <- exp(mean(log(phy$edge.length[phy$edge.length > 0] )) - 3 * sd(log(phy$edge.length[phy$edge.length > 0])))                          
phy$edge.length[phy$edge.length < time_min] <- time_min / 2

  
ratevec <- as.vector((tipdist/cophenetic(phy))[lower.tri(tipdist)])
                             
sigma_prior <- sqrt(mean(ratevec[ratevec > 0])) 

                             
data <- list(N = N,
             NI = NI,
             phi_mu = latlonrad[,1],
             theta_mu = latlonrad[,2],
             time = phy$edge.length,
             self = phy$edge[,2],
             ancestor = phy$edge[,1],
             sigma_prior = sigma_prior,
             kappa = latlon_mat[,3])

init <- list(sigma_raw = 1, loc_anc = latlon_inits_cartesian_noise)

write_stan_json(init, file.path(output_prefix, paste0(model_name, '_inits.json')))
write_stan_json(data, file.path(output_prefix, paste0(model_name, '_data.json')))

setwd(cmdstan_path())
system(paste0('make ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_command)
print(date())
system(sampling_command)

setwd(wd)
stanfit <- read_stan_csv(file.path(output_prefix, paste0(model_name, '_samples.txt')))
check_hmc_diagnostics(stanfit)
                             
summary(stanfit,pars='sigma_raw')
                             
                             
locations_cartesian_tensor <- extract(stanfit,pars='loc')[[1]]
locations_cartesian <- t(apply(locations_cartesian_tensor, 2, function(x) Directional:::vmf.mle(x)$mu))                         
locations_radians <- t(apply(locations_cartesian, 1, function(x) c(phi=atan2(x[2],x[1]),theta=atan2(sqrt(sum(x[1:2]^2)),x[3]))))
locations_degrees <- locations_radians / pi * 180
locations_degrees[,2] <- 90 - locations_degrees[,2]
colnames(locations_degrees) <- c('lon','lat')
rownames(locations_degrees) <- c(rownames(latlon_tips), paste0('ancestor_',(N+1):(N+NI)))                             
write.table(locations_degrees, file=file.path(output_prefix,'latlon_noise_tip_and_ancestors_north_america.txt'),sep='\t',quote=FALSE,row.names=FALSE)

library(leaflet)

#locdeg <- as.data.frame(locations_degrees)
locdeg <- as.data.frame(locations_degrees[1:N,])
#locdeg <- latlon_internal_guess
#colnames(locdeg) <- c('lon','lat')
#locdeg <- as.data.frame(latlon_mat)
#colnames(locdeg) <- c('lon','lat','kappa')
sp:::coordinates(locdeg) <- ~lon+lat
leaflet(locdeg) %>% addMarkers() %>% addTiles()
                             
sigma <- mean(extract(stanfit,pars='sigma_raw')[[1]])
