# manage dockerfiles

# Trim Galore
sudo docker build -t crpeters/trim_galore:0.6.7 -f dockerfiles/trim_galore/Dockerfile .
sudo docker push  crpeters/trim_galore:0.6.7

# Bedtools
sudo docker build -t crpeters/bedtools:latest -f dockerfiles/bedtools/Dockerfile .
sudo docker push  crpeters/bedtools:latest

# R-tidyverse-optparse
sudo docker build -t crpeters/r-tidyverse-optparse:4.1.2 -f dockerfiles/r-tidyverse-optparse/Dockerfile .
sudo docker push  crpeters/r-tidyverse-optparse:4.1.2
