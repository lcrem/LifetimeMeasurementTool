#!/bin/bash

firstRun=$1
lastRun=$2

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh



DIRECTORY=/data/PurityMonitor/Filling/


for runNum in `seq -w $firstRun $lastRun`;
do
#  echo $run
    run=$DIRECTORY/Run$runNum/

    if [ -f $run/RawAverages_ch4.root ]; then
	echo "Run is $runNum : calculating filtered avg for PrM 2"
	./findFilteredAverages $DIRECTORY $runNum 2
    fi
    
    if [  -f $run/RawAverages_ch1.root ]; then
	echo "Run is $runNum : calculating  avg for PrM 1"
	./findFilteredAverages $DIRECTORY $runNum 1
    fi

done
