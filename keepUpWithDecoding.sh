#!/bin/bash

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh

DIRECTORY=/data/PurityMonitor/Filling/


for run in $DIRECTORY/Run*
do
#  echo $run
  if [ ! -f $run/RootifiedFromBinary.root ]; then 

      if [ -f $run/Binary.bin ]; then 
	  ./decodeEvent $run 1
      fi 
  fi
  
done
