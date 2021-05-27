#!/bin/bash

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh

DIRECTORY=/data/PurityMonitor/Filling/


for run in $DIRECTORY/Run*
do
#  echo $run
    
    if [ -f $run/BadFlag ]; then
	continue;
    fi


  if [ -f $run/RawAverages_ch4.root ]; then 

      if [ ! -f $run/PrM2_filtAvg.root ]; then
	  runNumber=$(echo $run | sed 's/[^0-9]*//g')
	  echo "Run is $runNumber : filtered avg for PrM 2 not found!"
	  ./findFilteredAverages $DIRECTORY $runNumber 2
      fi
  fi

  if [ -f $run/RawAverages_ch1.root ]; then 
      
      if [ ! -f $run/PrM1_filtAvg.root ]; then
	  runNumber=$(echo $run | sed 's/[^0-9]*//g')
	  echo "Run is $runNumber : filtered avg for PrM 1 not found!"
	  ./findFilteredAverages $DIRECTORY $runNumber 1
      fi
  fi
  
done
