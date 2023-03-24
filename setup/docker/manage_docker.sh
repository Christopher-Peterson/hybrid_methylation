# manage dockerfiles
cd setup/docker

# Trim Galore
sudo docker build -t crpeters/trim_galore:0.6.7 -f trim_galore/Dockerfile .
sudo docker push  crpeters/trim_galore:0.6.7

# Bedtools
sudo docker build -t crpeters/bedtools:latest -f bedtools/Dockerfile .
sudo docker push  crpeters/bedtools:latest

# R-tidyverse-optparse
sudo docker build -t crpeters/r-tidyverse-optparse:4.2.1 -f r-tidyverse-optparse/Dockerfile .
sudo docker push  crpeters/r-tidyverse-optparse:4.2.1

# R-plotting
sudo docker build -t crpeters/r-plotting:4.2.1 -f r-plotting/Dockerfile .
sudo docker push  crpeters/r-plotting:4.2.1


# SNPsplit
sudo docker build -t crpeters/snpsplit:0.5.0 -f snpsplit/Dockerfile .
sudo docker push  crpeters/snpsplit:0.5.0

# ANGSD
sudo docker build -t crpeters/angsd:0.940-stable -f angsd/Dockerfile .
sudo docker push  crpeters/angsd:0.940-stable

# BatMeth2 
sudo docker build -t crpeters/batmeth2:latest -f batmeth2/Dockerfile .
sudo docker push  crpeters/batmeth2:latest
 
# MethHaplo

sudo docker build -t crpeters/methhaplo:latest -f methhaplo/Dockerfile .
sudo docker push  crpeters/methhaplo:latest

# InformME.jl and CpelAsm.jl (note; InformMe is not actually included in this right now)
sudo docker build -t crpeters/inform_me_jl:latest -f inform_me_jl/Dockerfile .
sudo docker push crpeters/inform_me_jl:latest

# R-methylkit
sudo docker build -t crpeters/r-methylkit:4.2.1 -f r-methylkit/Dockerfile . 
sudo docker push  crpeters/r-methylkit:4.2.1

# WhatsHap
sudo docker build -t crpeters/whatshap:1.7 -f whatshap/Dockerfile . 
sudo docker push  crpeters/whatshap:1.7


