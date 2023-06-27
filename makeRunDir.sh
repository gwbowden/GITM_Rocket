#!/bin/bash
#Script written by Melrose to set up GITM run directory in a more portable manner

echo "Setting up GITM run directory"
echo "Make sure you're in /short/n23"
echo "Usage: makeRunDir.sh <directory name> (default name: run"

TOPDIR=${PWD}

# Set the run directory name based on supplied argument
if [ -z $1 ]
    then 
     echo "No argument supplied, using default directory name"
     name=run
    else
     name=$1
fi

# Warn and exit if the run directory already exists.
if [ -d "$name" ]
    then
     if [ -L "$name" ]
         then
          echo "Directory already exists as a symbolic link - exiting"
          exit 1
     else
      echo "Directory exists - exiting"
      exit 1
     fi
fi

# Create run directory and symbolic links. 
# Copied more or less word for word from the rundir section of $GITM_HOME/Makefile 
mkdir -p ${TOPDIR}/$name/UA

cd ./$name
ln -s ${EMPIRICALIEDIR}/data ./EIE
ln -s ${BINDIR}/PostProcess.exe ./PostGITM.exe

cp ${GITM_HOME}/rungtim.pbs ${TOPDIR}/$name/ 

cd ${TOPDIR}/$name/UA
mkdir restartOUT data DataIn
ln -s restartOUT restartIN
ln -s ${BINDIR}/pGITM
ln -s ${GITM_HOME}/srcData/* DataIn/
rm -rf ./DataIn/CVS
ln -s ${GITM_HOME}/data/* DataIn/
rm -rf ./DataIn/CVS
cd ${TOPDIR}/$name
ln -s ${BINDIR}/GITM.exe ./
cp UA/DataIn/UAM.in ./
touch core 
chmod 444 core
ln -s UA/* ./
echo "Run directory set up"



