#!/bin/sh

# set -x

INPIT_N2W=./dataset/MERGE_DARTS_201903250500.nc
INPUT_MRT=./mrtV2_config.txt
INPUT_N2P=./nexrad_ref.png

PROGRAM_PATH=..
OUTPUT_DIR=./output

mkdir -p $OUTPUT_DIR/
rm -rf $OUTPUT_DIR/*



sorted_files=()

PREP_FILES () {
    STAGE_1_FILES=()

    # Iterate through files in the directory
    for file in "$OUTPUT_DIR"/*; do
        # Check if the current item is a file
        if [ -f "$file" ]; then
            # Add file to the array with its last update time
            STAGE_1_FILES+=("$file")
        fi
    done
    # # Sort the files based on their last update time
    # sorted_files=($(printf '%s\n' "${STAGE_1_FILES[@]}" | sort -n -k 2))
    # Sort the files based on creation time
    sorted_files=($(printf '%s\n' "${STAGE_1_FILES[@]}" | sort -k2 -n | cut -d' ' -f1))
}

NOWCAST () {
    ## Stage 1
    $PROGRAM_PATH/NowcastToWDSS2/NowcastToWDSS2 $INPIT_N2W $OUTPUT_DIR
}

MRTV2 () {
    # Stage 2
    # for file in "${sorted_files[@]}"
    for ((i=0; i<=30; i++))
    do
        file="PredictedReflectivity_${i}min_mrtV2_co-nfig00.nc"
        echo "$file"
        $PROGRAM_PATH/mrtV2/mrtV2 -c $INPUT_MRT $OUTPUT_DIR/$file
    done
}

NETCDF2PNG () {
    # Stage 3
    for ((i=0; i<=30; i++))
    do
        INPUT_PR_FILE="PredictedReflectivity_${i}min_mrtV2_co-nfig00.nc"
        OUTPUT_PR_FILE="PredictedReflectivity_${i}min_mrtV2_co-nfig00.png"
        $PROGRAM_PATH/netcdf2png/merged_netcdf2png -c $INPUT_N2P -q 235 -z 0,75 -o $OUTPUT_DIR/$OUTPUT_PR_FILE $OUTPUT_DIR/$INPUT_PR_FILE
    done
}



NOWCAST

PREP_FILES

MRTV2

# NETCDF2PNG










