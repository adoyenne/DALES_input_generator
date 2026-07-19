#!/bin/bash
set -e

# “rolling short-range forecast composite”
# First cycle: 0–27h (DALES spin-up OK)
# Other cycles: take 3–27 (25 hours)

echo $STARTDATE
echo $ENDDATE
echo $SRC_DIR_gribs
echo $PATH_GRIB_COMP

# Input configuration
START_DATE="${STARTDATE}00"
END_DATE="${ENDDATE}00"

# Derived values
SRC_DIR=$SRC_DIR_gribs
OUT_DIR=$PATH_GRIB_COMP
DATE_FORMAT="%Y%m%d%H%M"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

# Function to increment date in HH steps using Python
date_add_hours() {
    python3 -c "from datetime import datetime, timedelta; \
d=datetime.strptime('$1','%Y%m%d%H'); \
d+=timedelta(hours=$2); \
print(d.strftime('%Y%m%d%H'))"
}

# Generate list of run start times between START_DATE and END_DATE
echo "Generating forecast cycle list from $START_DATE to $END_DATE..."
DATES=()
current_date="$START_DATE"

while [[ "$current_date" < "$END_DATE" || "$current_date" == "$END_DATE" ]]; do
    DATES+=("$current_date")
    current_date=$(date_add_hours "$current_date" 24)
done

echo "Forecast cycles:"
printf ' - %s\n' "${DATES[@]}"

# Step 1: Unzip all needed zip files
for date in "${DATES[@]}"; do
    short=${date:0:10}
    echo "Unzipping $short..."
    archive=$(printf "$HA_archive_pattern" "$short")
    mkdir -p "$SRC_DIR/$short"
    case "$archive" in
        *.tar)
            tar -xf "$SRC_DIR/$archive" -C "$SRC_DIR/$short"
            ;;
        *.zip)
            unzip -q "$SRC_DIR/$archive" -d "$SRC_DIR/$short"
            ;;
        *)
            echo "ERROR: Unsupported archive format: $archive" >&2
            exit 1
            ;;
    esac
done

# Step 2: Process and copy files
echo "Copying and renaming files to $OUT_DIR..."

forecast_hour=0   # Starts from 0
missing_count=0   # track missing files

for i in "${!DATES[@]}"; do
    current_date="${DATES[$i]}"
    current_short=${current_date:0:10}
    input_dir="$SRC_DIR/$current_short"

    if [[ "$i" -eq 0 ]]; then
        # First cycle: take 0–27 (28 hours)
        start_h=0
        end_h=27
    else
        # Other cycles: take 3–27 (25 hours)
        start_h=3
        end_h=27
    fi

    for ((h=start_h; h<=end_h; h++)); do
        hhh=$(printf "%03d" "$h")
        input_file="$input_dir/$(printf "$HA_grib_file_pattern" "$current_short" "$hhh")"
        
        if [[ -f "$input_file" ]]; then
            new_hour=$(printf "%03d" $forecast_hour)
            outfile=$(printf "$HA_grib_file_pattern" "$START_DATE" "$new_hour")
            cp "$input_file" "$OUT_DIR/$outfile"
            forecast_hour=$((forecast_hour + 1))
        else
            echo "Missing file: $input_file" >&2
            missing_count=$((missing_count + 1))
        fi
    done

    # clean up unzipped directory AFTER processing each cycle
    rm -rf "$input_dir"

done

# Final report
expected_files=$((28 + ((${#DATES[@]} - 1) * 25)))
grib_find_pattern="${HA_grib_file_pattern//%s/*}"
actual_files=$(find "$OUT_DIR" -name "$grib_find_pattern" | wc -l)

echo "Copied $actual_files files to $OUT_DIR"

if [[ "$actual_files" -ne "$expected_files" ]]; then
    echo "Warning: expected $expected_files files, but found $actual_files." >&2
fi

# report missing files clearly
if [[ "$missing_count" -gt 0 ]]; then
    echo "ERROR: $missing_count files were missing!" >&2
    exit 1
fi

# sanity check on timeline consistency
if [[ $forecast_hour -ne $((actual_files)) ]]; then
    echo "Time indexing mismatch! Check for gaps." >&2
fi

echo "Total forecast hours: $((forecast_hour))"
echo "Done. Forecast window from $START_DATE to $END_DATE."