
# Installing IGV and IGVtools

We’ll be using [IGV](https://igv.org) and the associated command line
tools to sanity-check our alignments. First, download and unzip both IGV
(with Java) and the IGV + IGVtools bundle (without Java) to TACC.

``` bash
cdw 
VERSION=2.12.3
REMOTE=data.broadinstitute.org/igv/projects/downloads/2.12
# IGV w Java
wget https://${REMOTE}/IGV_Linux_${VERSION}_WithJava.zip 
# IGV tools
wget https://${REMOTE}/IGV_${VERSION}.zip

# Unzip and clean up
unzip IGV_Linux_${VERSION}_WithJava.zip
unzip IGV_${VERSION}.zip

# Move IGV tools to the directory that has  the bundled java
mv IGV_${VERSION}/igvtools* IGV_Linux_${VERSION}

# Cleanup
rm IGV*zip
rm -r IGV_${VERSION} # folder w/o java 
```

Next, we want to make sure that IGV doesn’t start writing data to $HOME
when we run it. We’ll do this by making a scratch directory for IGV and
then linking it to home.

``` bash
cdw
mkdir IGV_genomes
cdh
mkdir igv
cd igv
ln -s $WORK/IGV_files genomes
```

## Using IGV

We’re going to use the [TACC Visualization Portal](vis.tacc.utexas.edu/)
to use IGV. Go to the portal and start a DCV remote desktop session on
Lonestar 6 using the development queue. If DCV is unavailable, you can
use VNC instead. When the page changes, click the green “Connect” button
(once it appears), then sign into TACC on the DCV page.

On the virtual desktop, open a terminal, navigate to where you saved
IGV, and run `igv.sh`. If you want, you can right-click on the desktop
and create a launcher (a.k.a., a shortcut) for it to save time in the
future.
