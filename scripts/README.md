code to download and massage data, run models, and produce visuals for inclusion in the shiny app

aggregate_inat.r - downloads inaturalist observations, will include code to aggregate these observations into counts-per-county and produce GPS coordinates for the cetnroids of the counties

north-american-data-analysis.R - organize nextstrain phylogeny + gps data

run_alr.r - estimate gps coordinates of ancestral sequences in a phylogeny based on time-scaled edges and gps-located tips

run_gp_continuous.r - estimate value of a continuous variable at arbitrary points on the globe based on observations of that variable at other points on the globe

quick_alr.r - association location strings from NextStrain phylogenetic tips with GPS coordinates, then estimate ancestral locations of internal nodes in phylogeny with simple assumption - the centroid location of a polygon formed by all descendant tips

glmer_quick.r - simple statistical test for assocation between number of mutations and corresponding covariates

2/24/2021 : 4:30 PM
                                    Samuel Coleman uploaded the inital version of north-american-data-analysis
    The purpose of this file is to analyze data from NextStrain regarding the ongoing Covid Pandemic within North America. This script isolates the Covid sequences from NextStrain that have location data within their metadata.
    
