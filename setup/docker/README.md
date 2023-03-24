
# Docker/Singularity Setup

This project uses Docker (and Singularity, the TACC-equivalent version)
for installation and version control of bioniformatics packages. When
installed correctly, these are invisible to the end-user.

## Installation

1.  Setup your directories with `mkdir -p $WORK/singularity/bin`.
2.  Copy `docker/singularity_mask.sh` to `$WORK/singularity`.
3.  Edit your `.bashrc` to include `$WORK/singularity/bin` in your
    `$PATH`.
4.  Prepare and run the singularity installation script:

``` bash
cds hybrid_methylation/docker

mkdir logs
sbatch slurm/setup_singularity.slurm
```

This will install the singularity modules and configure them so that
they can be used as regular bash commands.

## Setting up the brms/stan docker image

The GAM components (the temporary ones, at least) use an image
Christopher originally created for a different project; it’s being set
up here for both local and TACC use.

``` bash
TAG=crpeters/edge-trait-meta:4.2.0 
sudo docker pull $TAG

# Configure the image:
# Set the working directory as your current directory (the project root)
LOCAL_DIR=$PWD 
CONF_DIR=$LOCAL_DIR/setup/docker/prefs
# If you'd like to access the git repo from within the image, mount your ssh directory
mkdir $CONF_DIR
# Otherwise, comment this line out
SSH_MOUNT="-v $HOME/.ssh:/home/rstudio/.ssh "
INTERNAL_NAME='hybrid_methylation'

NAME=r-brms
PORT=8784
# Run the container
sudo docker run -d -p 127.0.0.1:$PORT:8787 \
-v $LOCAL_DIR:/home/rstudio/${INTERNAL_NAME} \
-v $CONF_DIR:/home/rstudio/.config/rstudio \
$SSH_MOUNT \
-e DISABLE_AUTH=true -e USERID=$UID --name $NAME $TAG
```

Note that this will configure RStudio Server to run without
authentication. If you’re running this on a local machine, this is the
easiest option; however, running it on a resource that is accessible
from the web would be a security vulnerability. In this case, replace
`-e DISABLE_AUTH=true` with `-e PASSWORD=...` (where the dots are your
password).

### Using the brms Docker Image

In the console, run `sudo docker container start edge_trait_meta` (you
can skip this if you just initialized the container). Then, open a web
browser and go to [localhost:8784](localhost:8784) (Change the last
number to whatever \$PORT is if you adjust it).

When you’re done, run `sudo docker container stop edge_trait_meta` to
shut it down.

### If there are errors

There seems to be an issue with Docker-Hub where some of my newer
repositories not wanting to download on Singularity. This appears to be
some sort of privacy/authorization problem. Logging into my account
fixes it. If someone wants to download these images and gets any errors,
contact Christopher and he’ll give you a read-only access code.
