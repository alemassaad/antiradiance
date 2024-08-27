#!/bin/bash

LOG_FILE="timing_results.log"

#clear previous logs
>$LOG_FILE

echo "Processing legoBrick.obj with configuration lego_brick..." | tee -a $LOG_FILE
START_TIME=$(date +%s.%N)
./antiradiance "objects/legoBrick.obj" "lego_brick"
END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)
echo "Time taken for legoBrick.obj: $ELAPSED_TIME seconds" | tee -a $LOG_FILE
echo "-------------------------------------------" | tee -a $LOG_FILE

echo "Processing room_radiosity_small.obj with configuration small_room..." | tee -a $LOG_FILE
START_TIME=$(date +%s.%N)
./antiradiance "objects/room_radiosity_small.obj" "small_room"
END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)
echo "Time taken for room_radiosity_small.obj: $ELAPSED_TIME seconds" | tee -a $LOG_FILE
echo "-------------------------------------------" | tee -a $LOG_FILE

echo "Processing bear.obj with configuration bear..." | tee -a $LOG_FILE
START_TIME=$(date +%s.%N)
./antiradiance "objects/bear.obj" "bear"
END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)
echo "Time taken for bear.obj: $ELAPSED_TIME seconds" | tee -a $LOG_FILE
echo "-------------------------------------------" | tee -a $LOG_FILE

echo "All meshes processed. Timing results saved in $LOG_FILE." | tee -a $LOG_FILE
