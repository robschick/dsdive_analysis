#!/bin/bash

echo Updating github dsdive repos...

if [ -d "./dsdive" ];
then
  cd dsdive
  git pull
  cd ..
else
  git clone git@github.com:jmhewitt/dsdive.git
fi

echo Installing package into singularity image...

singularity exec rmovement.sif install2.r -l libs -r NULL ./dsdive

echo Finished updating dsdive installation.
