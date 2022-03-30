#!/bin/bash

# Rscript analysis/01-load-data.R

# cd data-generated;rm index-*
# cd ..

# # Rscript analysis/02-fit-geostat.R 'SYN' 'tweedie'
# Rscript analysis/02-fit-geostat.R 'SYN' 'binomial_gamma'
# # Rscript analysis/02-fit-geostat.R 'HBLL' 'tweedie'
# # Rscript analysis/02-fit-geostat.R 'HBLL' 'nbinom2'
# Rscript analysis/02-fit-geostat.R 'HBLL' 'binomial_gamma'

Rscript analysis/03-process-plot-geostat.R 'SYN' 'binomial_gamma'
Rscript analysis/03-process-plot-geostat.R 'HBLL' 'binomial_gamma'

Rscript analysis/04-rel-error-slopes.R 'SYN'
Rscript analysis/04-rel-error-slopes.R 'HBLL'

Rscript analysis/05-plot-all-survs.R 0
