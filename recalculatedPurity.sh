#!/bin/bash

firstRun=$1
lastRun=$2

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh



DIRECTORY=/data/PurityMonitor/Filling/


for runNum in `seq $firstRun $lastRun`;
do
#  echo $run
    run=$DIRECTORY/Run$runNum/

    if [ -f $run/PrM2_filtAvg.root ]; then
	echo "Run is $runNum : filtered avg for PrM 2 not found!"
	./calculateLifetimeWithDigitiser $DIRECTORY $runNum 2
    fi
    
    if [  -f $run/PrM1_filtAvg.root ]; then
	echo "Run is $runNum : filtered avg for PrM 1 not found!"
	./calculateLifetimeWithDigitiser $DIRECTORY $runNum 1
    fi

done
