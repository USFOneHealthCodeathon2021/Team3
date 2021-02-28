library(ape)
library(cmdstanr)
library(rstan)

wd <- getwd() # start in the Team 3 directory
output_prefix <- file.path(wd,'outputs/alr_gp')
model_dir <- file.path(wd, 'models')
model_name <- 'alr_gp'
dataset_name <- 'population'

sampling_command <- paste(paste0('./', model_name),
                paste0('data file=', file.path(output_prefix, paste0(model_name, '_', dataset_name, '_data.json'))),
                paste0('init=', file.path(output_prefix, paste0(model_name, '_', dataset_name, '_inits.json'))),
                'output',
                paste0('file=', file.path(output_prefix, paste0(model_name, '_', dataset_name, '_samples.txt'))),
                paste0('refresh=', 100),
                'method=sample',
                'adapt delta=0.8',
                'algorithm=hmc',
                #'stepsize=0.01',
                'engine=nuts',
                'max_depth=10',
                'num_warmup=500',
                'num_samples=200',
                sep=' ')

## import and deal with NextStrain dataset

phy_time <- collapse.singles(read.tree(file.path(wd,'datasets/raw/NextStrain/nextstrain_ncov_north-america_timetree.nwk')))
tipdat <- read.table(file.path(wd,'datasets/raw/NextStrain/nextstrain_ncov_north-america_metadata.tsv'),header=T, sep='\t', quote='')

tipdat$location_combined <- apply(tipdat[,c('Location','Admin.Division','Country')], 1, function(x) paste(x,collapse=', '))
#latlon_phy <- tmaptools:::geocode_OSM(unique(tipdat$location_combined))
#latlon_phy <- latlon_phy[!(is.na(latlon_phy$lat) | is.na(latlon_phy$lon)),]
#write.table(latlon_phy,file=file.path(wd,'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',quote=FALSE,row.names=FALSE)
latlon_phy <- read.table(file.path(getwd(),'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',header=TRUE)

latlon_phy$kappa <- apply(latlon_phy[,c('lon_min','lon_max','lat_min','lat_max')], 1, function(x) {
  r <- geosphere:::distGeo(c(x[[1]],x[[3]]),c(x[[2]],x[[4]]), a=1) / 2
  return(qnorm(0.99,0,1/sqrt(r)))
})

latlon_tips <- t(sapply(phy_time$tip.label, function(x) latlon_phy[latlon_phy$query==tipdat[tipdat$Strain==x,'location_combined'],c('lon','lat','kappa')]))                       
### temporary to scale down data
latlon_tips <- latlon_tips[log(as.numeric(latlon_tips[,'kappa'])) > 2.5,]
###

#phy_time <- drop.tip(phy_time,phy_time$tip.label[!phy_time$tip.label %in% rownames(latlon_tips)])

#write.tree(phy_time, file=file.path(output_prefix, 'nextstrain_ncov_north-america_timetree_filteredGPS.newick'))  
phy_time <- read.tree(file.path(output_prefix, 'nextstrain_ncov_north-america_timetree_filteredGPS.newick'))

latlon_mat <- matrix(as.numeric(latlon_tips[phy_time$tip.label,]),ncol=3)
latlonrad <- latlon_mat[,1:2] / 180 * pi
latlonrad[,2] <- (pi/2) - latlonrad[,2] ## geographic degrees go from -90 to +90 starting at equator, math starts at top and goes 0-pi
                        
latlon_internal_guess <- t(sapply(unique(phy_time$edge[,1]), function(x) icosa:::surfacecentroid(matrix(as.numeric(latlon_tips[phangorn:::Descendants(phy_time,x,'tips')[[1]],1:2]),ncol=2))))   
latlon_internal_guess_rad <- DescTools:::DegToRad(latlon_internal_guess)
latlon_internal_guess_rad[,2] <- (pi/2) - latlon_internal_guess_rad[,2]
                                  
latlon_inits_rad <- rbind(latlonrad, latlon_internal_guess_rad)
latlon_inits_cartesian <- t(apply(latlon_inits_rad, 1, function(x) c(sin(x[2])*cos(x[1]), sin(x[2])*sin(x[1]), cos(x[2])))) 
latlon_inits_cartesian_noise <- latlon_inits_cartesian + matrix(rnorm(nrow(latlon_inits_cartesian)*3,0,1e-4),ncol=3)
latlon_inits_cartesian_noise <- t(apply(latlon_inits_cartesian_noise, 1, function(x) x / sqrt(sum(x^2))))                                   
      
NT <- nrow(latlon_mat)
NI <- phy_time$Nnode

tipdist <- geosphere:::distm(latlon_mat[,1:2], fun = function(x,y) geosphere:::distGeo(x,y,a=1))

time_min <- exp(mean(log(phy_time$edge.length[phy_time$edge.length > 0] )) - 3 * sd(log(phy_time$edge.length[phy_time$edge.length > 0])))                          
phy_time$edge.length[phy_time$edge.length < time_min] <- time_min / 2

ratevec <- as.vector((tipdist/cophenetic(phy_time))[lower.tri(tipdist)])
                             
sigma_d_prior <- sqrt(mean(ratevec[ratevec > 0])) 

phy_div <- collapse.singles(read.tree(file.path(getwd(),'datasets/raw/NextStrain/nextstrain_ncov_north-america_tree.nwk')))
phy_div <- drop.tip(phy_div, phy_div$tip.label[!phy_div$tip.label %in% phy_time$tip.label])

###

### import and deal with covariate dataset

dat_co <- read.csv(file.path(wd, 'datasets/curated/census-app/Average_Household_Size_and_Population_Density_-_County.csv'),header=T)
dat_co$combined <- apply(dat_co,1,function(x) paste(x[6],x[7],sep=', '))
#latlon_co <- tmaptools:::geocode_OSM(dat_co$combined)
#write.table(latlon_co, file=file.path(wd,'datasets','curated','census-app','latlon_Average_Household_Size_and_Population_Density_-_County.txt'),sep='\t',quote=FALSE,row.names=FALSE)
latlon_co <- read.table(file.path(wd,'datasets','curated','census-app','latlon_Average_Household_Size_and_Population_Density_-_County.txt'),sep='\t',header=TRUE)
latlon_co <- latlon_co[!duplicated(latlon_co[,-1]),]
dat_co <- dat_co[dat_co$combined %in% latlon_co$query,]
latlon_co <- latlon_co[latlon_co$query %in% dat_co$combined,]

NY <- nrow(latlon_co)

y_log <- log(dat_co$B01001_calc_PopDensity)                     
scale_gp_prior <- sd(y_log)    

latlon_co_mat <- matrix(as.numeric(as.matrix(latlon_co[,c(3,2)])),ncol=2)
latlonrad_co <- latlon_co_mat[,1:2] / 180 * pi
latlonrad_co[,2] <- (pi/2) - latlonrad_co[,2] ## geographic degrees go from -90 to +90 starting at equator, math starts at top and goes 0-pi      
latlon_co_cartesian <- t(apply(latlonrad_co, 1, function(x) c(sin(x[2])*cos(x[1]), sin(x[2])*sin(x[1]), cos(x[2])))) 

###

### full dataset quantities

dists <- geosphere:::distm(rbind(latlon_phy[,c(3,2)],latlon_co[,c(3,2)]), fun = function(x,y) geosphere:::distGeo(x,y,a=1))

###
                             
data <- list(NT = NT,
             NI = NI,
             NY = NY,
             loc_mu = t(latlon_inits_cartesian[1:NT,]),
             kappa = latlon_mat[,3],
             time = phy_time$edge.length,
             self = phy_time$edge[,2],
             ancestor = phy_time$edge[,1],
             sigma_d_prior = sigma_d_prior,
             loc_y = t(latlon_co_cartesian),
             y_obs = y_log,
             mut = phy_div$edge.length,
             rho_prior = mean(dists[lower.tri(dists)]),
             scale_gp_prior = scale_gp_prior)

init <- list(rho_raw = 0.1, locvec = t(latlon_inits_cartesian_noise))

write_stan_json(init, file.path(output_prefix, paste0(model_name, '_', dataset_name, '_inits.json')))
write_stan_json(data, file.path(output_prefix, paste0(model_name, '_', dataset_name, '_data.json')))

setwd(cmdstan_path())
system(paste0('make ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_command)
print(date())
system(sampling_command)

setwd(wd)
stanfit <- read_stan_csv(file.path(output_prefix, paste0(model_name, '_', dataset_name, '_samples.txt')))
check_hmc_diagnostics(stanfit)
                             
summary(stanfit,pars='sigma_d_raw')
                             
                             
locations_cartesian_tensor <- extract(stanfit,pars='loc')[[1]]
locations_cartesian <- t(apply(locations_cartesian_tensor, 3, function(x) Directional:::vmf.mle(x)$mu))                         
locations_radians <- t(apply(locations_cartesian, 1, function(x) c(phi=atan2(x[2],x[1]),theta=atan2(sqrt(sum(x[1:2]^2)),x[3]))))
locations_degrees <- locations_radians / pi * 180
locations_degrees[,2] <- 90 - locations_degrees[,2]
colnames(locations_degrees) <- c('lon','lat')
rownames(locations_degrees) <- c(rownames(latlon_tips), paste0('ancestor_',(NT+1):(NT+NI)))                             
write.table(locations_degrees, file=file.path(output_prefix, paste0(model_name, '_', dataset_name, 'latlon_noise_tip_and_ancestors_north_america.txt')),sep='\t',quote=FALSE,row.names=FALSE)

library(leaflet)

#locdeg <- as.data.frame(locations_degrees)
locdeg <- as.data.frame(locations_degrees[1:NT,])
#locdeg <- latlon_internal_guess
#colnames(locdeg) <- c('lon','lat')
#locdeg <- as.data.frame(latlon_mat)
#colnames(locdeg) <- c('lon','lat','kappa')
sp:::coordinates(locdeg) <- ~lon+lat
leaflet(locdeg) %>% addMarkers() %>% addTiles()
                             
sigma_d <- mean(extract(stanfit,pars='sigma_d_raw')[[1]])
