# Project Name: Viral SpaceTime
---
**Team Leaders**: Ryan McMinds, Jesper Madsen

**Team Members**:  Ann Mathew, Sofia	Bhatia, Samuel Coleman, Omkar	Dokur, Janelle Donglasan, Jan Dahrendorff

**GVN/USF mentors**: Dr. Alash'le G. Abimiku



## Objectives

Linking geospatial datasets with viral evolution rates

## Methods and Implementation

1st pass: 

  - download NextStrain US subset phylogenies, including time-scaled and divergence-scaled branch lengths, and geographic locations of tips.
  - impute mean geographic location of edges. 
  -   first impute ancestral nodes as a point on the shortest great-circle arc between the direct descendant nodes/tips. start at tips and iterate toward root of tree. The point for each ancestor corresponds to the location proportional to time-scaled branch lengths (e.g. if one descendant edge has length 2, and the other has length 1, the point is 1/3 the distance from descendant 1 to descendant 2.
  -   then impute mean location for edges as midpoint between the ancestor and descendant
  -   input imputed locations/distances into gaussian process or inla model

- inla/gp model conditioned on geographic dataset, with imputed values at locations corresponding to the estimated location of phylogenetic edges


![Flowchart](https://github.com/USFOneHealthCodeathon2021/Team3/blob/main/images/flowchart-v1.pdf)


## Results 



## Running our App

The following packages are required to run Viral SpaceTime:

RStudio

R-packages:
 * library(shiny)
 * library(leaflet)
 * library(leaflet.extras)
 * library(dplyr)
 * library(DT)

Download and open the "app.R" file in RStudio. In the upper-right corner of RStudio, press the "Run App" button.

![RStudio screenshot](https://github.com/USFOneHealthCodeathon2021/Team3/blob/main/images/interactive_map_1.png)

![RStudio screenshot](https://github.com/USFOneHealthCodeathon2021/Team3/blob/main/images/interactive_map_2.png)

![RStudio screenshot](https://github.com/USFOneHealthCodeathon2021/Team3/blob/main/images/phy_data.png)

![RStudio screenshot](https://github.com/USFOneHealthCodeathon2021/Team3/blob/main/images/test.png)
