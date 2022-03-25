#!/bin/bash

Rscript analysis/01-load-data.R

cd data-generated;rm index-*
# Rscript analysis/02-fit-geostat.R 'SYN' 'tweedie'
Rscript analysis/02-fit-geostat.R 'SYN' 'binomial_gamma'
# Rscript analysis/02-fit-geostat.R 'HBLL' 'tweedie'
# Rscript analysis/02-fit-geostat.R 'HBLL' 'nbinom2'
Rscript analysis/02-fit-geostat.R 'HBLL' 'binomial_gamma'

