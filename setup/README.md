
# Installation and Setup

These sub-directories contain instructions for various setup and data
preparation tasks.

- [The Singularity images needed for the workflow](docker) (Do this
  first!)
- [The WGBS data and genome](wgbs_setup)
- [Integrative Genomics Viewer on TACC](igv)
- [Additional setup for contributers to this repository](dev)

## Julia Setup

Given the difficulty in getting Julia to work wtih a docker image, I’m
just going to try running it natively. This project requires Julia 1.3.1
to run CpelAsm/InformME.jl

### Install Julia

``` bash
cdw
wget https://julialang-s3.julialang.org/bin/linux/x64/1.3/julia-1.3.0-linux-x86_64.tar.gz
tar -xvf julia-1.3.0-linux-x86_64.tar.gz
rm julia-1.3.0-linux-x86_64.tar.gz
cd julia-1.3.0
# Move .julia from home to work
mv $HOME/.julia .
cdh 
ln -s $WORK/julia-1.3.0/.julia .julia
# Add t path
ln -s $WORK/julia-1.3.0/bin/julia $HOME/bin/julia-1.3.0
```

### Install CpemASM

``` bash

julia-1.3.0
] # Opens Pkg interface
add https://github.com/jordiabante/CpelAsm.jl.git
test CpelAsm
# backspace to get out of pkg>
using CpelAsm
```

### Install InformME.jl

``` bash
julia-1.3.0
] # open Pkg
# Install dependency
add https://github.com/timholy/QuadDIRECT.jl.git
# Install package
add https://github.com/GarrettJenkinson/InformMe.jl
test InformMe
# backspace to get out of pkg>
using InformMe
```
