library(ape)

phy <- read.tree('/Users/Ryan/codeathon2021/Team3/datasets/raw/NextStrain/nextstrain_ncov_north-america_timetree.nwk')

dat <- read.table('/Users/Ryan/codeathon2021/Team3/datasets/raw/NextStrain/nextstrain_ncov_north-america_metadata.tsv',header=T, sep='\t', quote='')

dat <- dat[dat$Location != '',]

phy <- drop.tip(phy,phy$tip.label[!phy$tip.label %in% dat$Strain])

library('tmaptools')
#latlon <- geocode_OSM(unique(dat$Location))
#write.table(latlon,file=file.path(getwd(),'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',quote=FALSE,row.names=FALSE)
read.table(file.path(getwd(),'datasets','curated','NextStrain','latlon_north_america.txt'),sep='\t',header=TRUE)

library(geosphere)
gcIntermediate(latlon[1,c(3,2)],latlon[2,c(3,2)],1)

latlon_tips <- t(sapply(phy$tip.label,function(x) latlon[latlon$query==dat[dat$Strain==x,'Location'],c('lon','lat')]))

                        library(phangorn)
                        library(icosa)
latlon_internal <- t(sapply(unique(phy$edge[,1]), function(x) surfacecentroid(matrix(as.numeric(latlon_tips[Descendants(phy,x,'tips')[[1]],]),ncol=2)))  )   
                            
alllatlon <- rbind(latlon_tips,latlon_internal)
                            
write.table(alllatlon,file=file.path(getwd(),'outputs','quick_alr','latlon_north_america_phy.txt'),sep='\t',quote=FALSE,row.names=FALSE)
write.tree(phy,file=file.path(getwd(),'outputs','quick_alr','phy_filtered.newick'))
