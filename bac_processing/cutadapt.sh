#!/bin/bash

# $1 = path to fastq directory
# $2 = forward primers
# $3 = reverse primers

# Changing input field separators to new line (\n)
IFS=$'\n'

# Create a variable for the current path
current=$PWD

# Error message for missing positional arguments when running this script
if [ $# -ne 3 ] ; then
    echo "Usage:"
    echo "$0 takes three parameters"
    echo '$1 = path to fastq directory'
    echo '$2 = forward primers'
    echo '$3 = reverse primers'
    exit 1
else
    echo 'cutadapt batch run'
fi

# Adding a path to fastq directory from the positional argument $1 
fastqpath=$1
# Adding fwd primer sequences from the positional argument $2
fwdprm=$2
# Adding rvs primer sequences from the positional argument $3
rvsprm=$3

# Extracting fwd fastq file names
fwd=($(ls -1v $1 | grep 'R1'))

# Creating a new directory for clipped fastq files
mkdir $1/../Clipped 

# Move to fastq directory
cd $1

# Running cutadapt using forloop for each pair-end fastq files. 
for rvsfastq in "${fwd[@]}"
do
	rvsname=$(echo "$rvsfastq" | sed -e "s/R1/R2/")
	echo "$rvsfastq" >> $current/cutadapt_summary.txt
	cutadapt -g "$fwdprm" -G "$rvsprm" --discard-untrimmed -m 10 -o ./../Clipped/"$rvsfastq" -p ./../Clipped/"$rvsname" "$rvsfastq" "$rvsname" >> $current/cutadapt_summary.txt
done