#!/bin/bash

echo Downloading singularity image...

curl -O https://research-singularity-registry.oit.duke.edu/jmh182/rmovement.sif


echo Making libs directory...

if [ ! -d "./libs" ];
then
  mkdir libs
fi


echo Updating github repos...

if [ -d "./composr" ];
then
  cd composr
  git pull
  cd ..
else
  git clone git@github.com:jmhewitt/composr.git
fi

if [ -d "./dsdive" ];
then
  cd dsdive
  git pull
  cd ..
else
  git clone git@github.com:jmhewitt/dsdive.git
fi


echo Installing packages into singularity image...

singularity exec rmovement.sif install2.r -l libs Rmpi
singularity exec rmovement.sif install2.r -l libs -r NULL ./snow_custom_error_dump
singularity exec rmovement.sif install2.r -l libs pracma
singularity exec rmovement.sif install2.r -l libs yaml
singularity exec rmovement.sif install2.r -l libs expm
singularity exec rmovement.sif install2.r -l libs bigmemory
singularity exec rmovement.sif install2.r -l libs Rdsm
singularity exec rmovement.sif install2.r -l libs -r NULL ./composr
singularity exec rmovement.sif install2.r -l libs -r NULL ./dsdive


echo Updating RMPI snow profile...

sed -i 's/\(library([A-Za-z0-9]*\)\()\)/\1, lib.loc = c("singularity\/libs", ".", .libPaths())\2/g' libs/snow/RMPISNOWprofile

echo Finished setting up singularity image
