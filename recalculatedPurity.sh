#!/bin/bash

firstRun=$1
lastRun=$2
noiseRun=$3

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh



DIRECTORY=/data/PurityMonitor/Filling/


for runNum in `seq -w $firstRun $lastRun`;
do
#  echo $run
    run=$DIRECTORY/Run$runNum/


    if [ -f $run/BadFlag ]; then
	continue;
    fi

    if [ -f $run/Calibration ]; then
	continue;
    fi



    if [ -f $run/PrM2_filtAvg.root ]; then
	echo "Run is $runNum : filtered avg for PrM 2 not found!"
	./calculateLifetimeWithDigitiser $DIRECTORY $runNum 2 $noiseRun
    fi
    
    if [  -f $run/PrM1_filtAvg.root ]; then
	echo "Run is $runNum : filtered avg for PrM 1 not found!"
	./calculateLifetimeWithDigitiser $DIRECTORY $runNum 1 $noiseRun
    fi

done
