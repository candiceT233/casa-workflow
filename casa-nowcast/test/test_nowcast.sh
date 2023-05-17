#!/bin/sh

# set -x

INPUT=./mrtV2_config.txt
INPUT_PNG=./nexrad_ref.png
PROGRAM_PATH=..
OUTPUT_DIR=./output

mkdir -p $OUTPUT_DIR/
rm -rf $OUTPUT_DIR/*

## Stage 1
$PROGRAM_PATH/NowcastToWDSS2/NowcastToWDSS2 $INPUT $OUTPUT_DIR

STAGE_1_FILES=()

# # Find all files in the directory and add them to the array in order of creation time
# while IFS= read -r -d '' file; do
#     STAGE_1_FILES+=("$file")
# done < <(find "$OUTPUT_DIR" -maxdepth 1 -type f -printf "%T@ %p\0" | sort -z -n | cut -d' ' -f2-)


# Iterate through files in the directory
for file in "$OUTPUT_DIR"/*; do
    # Check if the current item is a file
    if [ -f "$file" ]; then
        # Add file to the array with its last update time
        STAGE_1_FILES+=("$file")
    fi
done
# Sort the files based on their last update time
sorted_files=($(printf '%s\n' "${STAGE_1_FILES[@]}" | sort -n -k 2))
# # Sort the files based on creation time
# sorted_files=($(printf '%s\n' "${STAGE_1_FILES[@]}" | sort -k2 -n | cut -d' ' -f1))




# # Stage 2
# # for file in "${sorted_files[@]}"
# for ((i=0; i<=30; i++))
# do
#     file="PredictedReflectivity_${i}min_mrtV2_co-nfig00.nc"
#     echo "$file"
#     $PROGRAM_PATH/mrtV2/mrtV2 -c $INPUT $OUTPUT_DIR/$file
# done



# # Stage 3
# for ((i=0; i<=30; i++))
# do
#     INPUT_PR_FILE="PredictedReflectivity_${i}min_mrtV2_co-nfig00.nc"
#     OUTPUT_PR_FILE="PredictedReflectivity_${i}min_mrtV2_co-nfig00.png"
#     $PROGRAM_PATH/netcdf2png/merged_netcdf2png -c $INPUT_PNG -q 235 -z 0,75 -o $OUTPUT_DIR/$OUTPUT_PR_FILE $OUTPUT_DIR/$INPUT_PR_FILE
# done



