#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/inputFileList.txt /path/to/outputFileName.root"
    exit 1
fi

# Get the input arguments
inputFileList=$1
outputFileName=$2

# Run the ROOT macro with the provided arguments
root -l -q "runPairDst.C(\"${inputFileList}\", \"${outputFileName}\")"
