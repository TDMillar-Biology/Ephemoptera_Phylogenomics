#!/bin/bash


while read line
do echo 'testing' + (awk -F "\"*,\"*" '{print $1}')

done < ./Plate_3_key.csv

