#!/bin/bash
#
# Script to concatenate daily ERA5 tracking files into monthly files
# Created: October 27, 2025
#
# Usage: bash concat_era5_daily_to_monthly.sh
#

set -e  # Exit on error

# Activate environment with ncrcat
echo "Activating e3sm_unified environment..."
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# Define directories
IN_DIR="/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/era5_tracking"
OUT_DIR="/pscratch/sd/w/wcmca1/hackathon/era5"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

echo "Input directory: ${IN_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Define file prefixes
declare -a PREFIXES=("AR_tracks_era5_ll025sc" "TC_tracks_era5_ll025sc" "ETC_test_tracks_era5_ll025sc")

# Get unique year-month combinations from files
# Extract YYYYMM from filenames like: AR_tracks_era5_ll025sc.2019080100_2019080123.nc
echo "Scanning for available year-months..."
YEAR_MONTHS=$(ls ${IN_DIR}/AR_tracks_era5_ll025sc.*.nc 2>/dev/null | \
              sed -E 's/.*\.([0-9]{6})[0-9]{4}_[0-9]{10}\.nc/\1/' | \
              sort -u)

if [ -z "$YEAR_MONTHS" ]; then
    echo "ERROR: No input files found in ${IN_DIR}"
    exit 1
fi

echo "Found year-months: $(echo $YEAR_MONTHS | tr '\n' ' ')"
echo ""

# Process each prefix (AR, TC, ETC)
for PREFIX in "${PREFIXES[@]}"; do
    echo "=========================================="
    echo "Processing: ${PREFIX}"
    echo "=========================================="
    
    # Process each year-month
    for YYYYMM in ${YEAR_MONTHS}; do
        echo ""
        echo "Processing ${PREFIX} for ${YYYYMM}..."
        
        # Find all daily files for this month
        INPUT_FILES="${IN_DIR}/${PREFIX}.${YYYYMM}*.nc"
        OUTPUT_FILE="${OUT_DIR}/${PREFIX}.${YYYYMM}.nc"
        
        # Count files
        N_FILES=$(ls ${INPUT_FILES} 2>/dev/null | wc -l)
        
        if [ ${N_FILES} -eq 0 ]; then
            echo "  WARNING: No files found for ${PREFIX} ${YYYYMM}, skipping..."
            continue
        fi
        
        echo "  Found ${N_FILES} daily files"
        echo "  Output: ${OUTPUT_FILE}"
        
        # Check if output already exists
        if [ -f "${OUTPUT_FILE}" ]; then
            echo "  WARNING: Output file already exists, skipping..."
            echo "           Delete ${OUTPUT_FILE} to regenerate"
            continue
        fi
        
        # Two-step process:
        # 1. Convert first file to have unlimited time dimension (using ncks)
        # 2. Concatenate all files using ncrcat
        
        # Get list of files as array
        FILES=(${INPUT_FILES})
        FIRST_FILE="${FILES[0]}"
        TEMP_DIR="${OUT_DIR}/tmp_${PREFIX}_${YYYYMM}"
        mkdir -p "${TEMP_DIR}"
        
        echo "  Step 1: Converting ${N_FILES} files to have unlimited time dimension..."
        
        # Convert each file to have unlimited time dimension
        CONVERTED_FILES=""
        for ((i=0; i<${N_FILES}; i++)); do
            INPUT_FILE="${FILES[$i]}"
            BASENAME=$(basename "${INPUT_FILE}")
            TEMP_FILE="${TEMP_DIR}/${BASENAME}"
            
            # Show progress every 10 files
            if [ $((i % 10)) -eq 0 ] || [ $i -eq $((N_FILES - 1)) ]; then
                echo "    Converting file $((i+1))/${N_FILES}..."
            fi
            
            # Use ncks to make time dimension unlimited
            ncks -O --mk_rec_dmn time "${INPUT_FILE}" "${TEMP_FILE}" 2>/dev/null
            CONVERTED_FILES="${CONVERTED_FILES} ${TEMP_FILE}"
        done
        
        echo "  Step 2: Concatenating files..."
        # Now concatenate the converted files
        ncrcat -O -h ${CONVERTED_FILES} "${OUTPUT_FILE}"
        
        # Clean up temporary files
        rm -rf "${TEMP_DIR}"
        
        # Check success
        if [ $? -eq 0 ]; then
            # Get file size for verification
            SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)
            echo "  SUCCESS: Created ${OUTPUT_FILE} (${SIZE})"
        else
            echo "  ERROR: ncrcat failed for ${PREFIX} ${YYYYMM}"
            exit 1
        fi
    done
    
    echo ""
done

echo ""
echo "=========================================="
echo "Concatenation complete!"
echo "=========================================="
echo ""
echo "Output files created in: ${OUT_DIR}"
echo ""
echo "File summary:"
ls -lh "${OUT_DIR}"/*.nc 2>/dev/null || echo "No output files found"
echo ""
echo "To use these files, update your Python script:"
echo "  in_dir = \"${OUT_DIR}\""
echo "  basename_ar = \"AR_tracks_era5_ll025sc.*.nc\""
echo "  basename_tc = \"TC_tracks_era5_ll025sc.*.nc\""
echo "  basename_etc = \"ETC_test_tracks_era5_ll025sc.*.nc\""
echo ""
