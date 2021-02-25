read.table('quick_alr_gps')
read_tree(phy_filtered)
read_tree(phy_divergence)


glmer(phy_divergence$edge.length ~ log(phy$edge.length) + covariate, link='poisson')
