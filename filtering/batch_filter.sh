#!/bin/bash

pdb_directory="/Users/Adam/Desktop/ya/kv/1.3/late-nov-run/de-novo/outputs-KV13-NovR-dn/"
filtered_files_directory="$pdb_directory/filtered"
log_file="$filtered_files_directory/log.txt"

tds_script="/Users/Adam/adiachem/metrics/terminal_distance_score.py"
sesd_script="/Users/Adam/adiachem/metrics/sesd.py"
wsc_script="/Users/Adam/adiachem/metrics/wavelet_substructure_count.py"

# General Settings
binder_chain="A"
atoms="all" # Atoms to consider when running SESD and WSSC

# N and C Terminus Scoring Settings
reference_atoms="B,122,O:B,232,O:B,342,O:B,452,O"

# Smallest Enclosing Sphere Diameter Settings
n_samples=30
z_score_threshold=1

# Wavelet Secondary Structure Counting Settings
resolution=128
prominence=5
atoms="all"

# Filtering Thresholds
n_terminus_score_threshold=1.0
c_terminus_score_threshold=1.0
smallest_enclosing_sphere_diameter_threshold=45.0
num_secondary_structures_threshold=2

mkdir -p "$filtered_files_directory"
csv_file="$filtered_files_directory/output.csv"
header="filename,num_secondary_structures,smallest_enclosing_sphere_diameter,n_terminus_score,c_terminus_score"

# Function to append log entries
append_log() {
    log_entry="$1"
    log_entry="${log_entry%%$'\r'}"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $log_entry" >> "$log_file"
}

# Function for Terminus Scoring
perform_terminus_scoring() {
    start_time=$SECONDS
    total_iterations=$(find "$pdb_directory" -maxdepth 1 -type f -name "*.pdb" | wc -l)
    current_iteration=0

    for pdb_file in "$pdb_directory"/*.pdb; do
        filename=$(basename -- "$pdb_file")
        filename_no_ext="${filename%.*}"
        output=$(python3 "$tds_script" -i "$pdb_file" -c $reference_atoms -b $binder_chain)
        n_terminus_score=$(echo "$output" | awk -F 'N Term Score = |, C Term Score = ' '/N Term Score = [0-9.]+/{print $2}')
        c_terminus_score=$(echo "$output" | awk -F ', C Term Score = ' '{print $2}')

        # Check if the row exists in the CSV file
        if grep -q "^$filename_no_ext," "$csv_file"; then
            # If the row exists, update the existing row
            awk -v filename="$filename_no_ext" -v n_score="$n_terminus_score" -v c_score="$c_terminus_score" \
                -F ',' '$1 == filename { $4 = n_score; $5 = c_score; print; next } 1' OFS=',' "$csv_file" > "$csv_file.tmp"
            mv "$csv_file.tmp" "$csv_file"
        else
            # If the row doesn't exist, append a new row
            echo "$filename_no_ext,,,,$n_terminus_score,$c_terminus_score" >> "$csv_file"
        fi

        ((current_iteration++))
        progress=$((current_iteration * 100 / total_iterations))
        elapsed_time=$((SECONDS - start_time))
        ncts_progress="Scoring Terminals: [$current_iteration/$total_iterations] $progress% | Elapsed Time: $elapsed_time seconds\r"
        echo -ne "\033[0K\r$ncts_progress"

        # Append the last call of echo to the log
        if [ "$current_iteration" -eq "$total_iterations" ]; then
            append_log "$ncts_progress"
        fi
    done
    echo
}

# Function for Substructure Counting
perform_substructure_counting() {
    start_time=$SECONDS
    total_iterations=$(find "$pdb_directory" -maxdepth 1 -type f -name "*.pdb" | wc -l)
    current_iteration=0

    for pdb_file in "$pdb_directory"/*.pdb; do
        filename=$(basename -- "$pdb_file")
        filename_no_ext="${filename%.*}"
        output=$(python3 "$wsc_script" -i "$pdb_file" -r "$resolution" -p "$prominence" -a "$atoms")
        num_substructures=$(echo "$output" | awk -F 'n = |,' '/n = [0-9]+/{print $2}')

        # Check if the row exists in the CSV file
        if grep -q "^$filename_no_ext," "$csv_file"; then
            # If the row exists, update the existing row
            awk -v filename="$filename_no_ext" -v substructures="$num_substructures" \
                -F ',' '$1 == filename { $2 = substructures; print; next } 1' OFS=',' "$csv_file" > "$csv_file.tmp"
            mv "$csv_file.tmp" "$csv_file"
        else
            # If the row doesn't exist, append a new row
            echo "$filename_no_ext,$num_substructures,,,," >> "$csv_file"
        fi

        ((current_iteration++))
        progress=$((current_iteration * 100 / total_iterations))
        elapsed_time=$((SECONDS - start_time))
        wssc_progress="Counting Substructures: [$current_iteration/$total_iterations] $progress% | Elapsed Time: $elapsed_time seconds\r"
        echo -ne "$wssc_progress"

        # Append the last call of echo to the log
        if [ "$current_iteration" -eq "$total_iterations" ]; then
            append_log "\033[0K\r$wssc_progress"
        fi
    done
    echo
}

# Function for SESD Calculation
perform_sesd_calculation() {
    start_time=$SECONDS
    total_iterations=$(find "$pdb_directory" -maxdepth 1 -type f -name "*.pdb" | wc -l)
    current_iteration=0

    for pdb_file in "$pdb_directory"/*.pdb; do
        filename=$(basename -- "$pdb_file")
        filename_no_ext="${filename%.*}"

        # Run the SESD calculation script
        sesd_output=$(python3 "$sesd_script" -i "$pdb_file" -c $binder_chain -n $n_samples -z $z_score_threshold -a "$atoms")
        sesd_value=$(echo "$sesd_output" | awk -F 'SESD = |,' '{printf "%.3f", $2}')

        # Check if the row exists in the CSV file
        if grep -q "^$filename_no_ext," "$csv_file"; then
            # If the row exists, update the existing row
            awk -v filename="$filename_no_ext" -v sesd="$sesd_value" \
                -F ',' '$1 == filename { $3 = sesd; print; next } 1' OFS=',' "$csv_file" > "$csv_file.tmp"
            mv "$csv_file.tmp" "$csv_file"
        else
            # If the row doesn't exist, append a new row
            echo "$filename_no_ext,$num_substructures,$sesd_value,," >> "$csv_file"
        fi

        ((current_iteration++))
        progress=$((current_iteration * 100 / total_iterations))
        elapsed_time=$((SECONDS - start_time))
        sesd_progress="Measuring Diameters: [$current_iteration/$total_iterations] $progress% | Elapsed Time: $elapsed_time seconds\r"
        echo -ne "\033[0K\r$sesd_progress"

        # Append the last call of echo to the log
        if [ "$current_iteration" -eq "$total_iterations" ]; then
            append_log "$sesd_progress"
        fi
    done
    echo
}

# Function to filter and move PDB files based on score thresholds
filter_files() {
    start_time=$SECONDS
    total_iterations=$(find "$pdb_directory" -maxdepth 1 -type f -name "*.pdb" | wc -l)
    current_iteration=0

    for pdb_file in "$pdb_directory"/*.pdb; do
        filename=$(basename -- "$pdb_file")
        filename_no_ext="${filename%.*}"

        # Get terminus scores from CSV file
        n_terminus_score=$(awk -F',' -v filename="$filename_no_ext" '$1 == filename { print $4 }' "$csv_file")
        c_terminus_score=$(awk -F',' -v filename="$filename_no_ext" '$1 == filename { print $5 }' "$csv_file")

        # Get SESD value
        sesd_value=$(awk -F',' -v filename="$filename_no_ext" '$1 == filename { print $3 }' "$csv_file")

        # Get number of substructures
        num_substructures=$(awk -F',' -v filename="$filename_no_ext" '$1 == filename { print $2 }' "$csv_file")

        # Debug output
        # echo "Debug: filename=$filename, n_terminus_score=$n_terminus_score, c_terminus_score=$c_terminus_score, sesd_value=$sesd_value, num_substructures=$num_substructures"

        if [ -n "$n_terminus_score" ] && [ -n "$c_terminus_score" ] && [ -n "$sesd_value" ] && [ -n "$num_substructures" ]; then
            if [ "$(echo "$n_terminus_score > $n_terminus_score_threshold" | bc -l)" -eq 1 ] &&
               [ "$(echo "$c_terminus_score > $c_terminus_score_threshold" | bc -l)" -eq 1 ] &&
               [ "$(awk -v s="$sesd_value" -v t="$smallest_enclosing_sphere_diameter_threshold" 'BEGIN { print (s < t) }')" -eq 1 ] &&
               [ "$(awk -v num="$num_substructures" -v t="$num_secondary_structures_threshold" 'BEGIN { print (num > t) }')" -eq 1 ]; then

                # Copy PDB file
                cp "$pdb_file" "$filtered_files_directory"

                # Copy PNG file with the same filename if it exists
                png_file="$pdb_directory/$filename_no_ext.png"
                cp "$png_file" "$filtered_files_directory"
            fi
        else
            echo
            echo "One or more variables are not set."
        fi

        ((current_iteration++))
        progress=$((current_iteration * 100 / total_iterations))
        elapsed_time=$((SECONDS - start_time))
        filtering_progress="Filtering: [$current_iteration/$total_iterations] $progress% | Elapsed Time: $elapsed_time seconds"
        echo -ne "\033[0K\r$filtering_progress"
        # Append the last call of echo to the log
        if [ "$current_iteration" -eq "$total_iterations" ]; then
            append_log "$filtering_progress"
        fi
    done
    echo

    total_elapsed_time=$((SECONDS - start_time))
    filtered_files_count=$(find "$filtered_files_directory" -maxdepth 1 -type f -name "*.pdb" | wc -l)
    echo "Filtering Results: $filtered_files_count / $total_iterations ($((filtered_files_count * 100 / total_iterations))%)"
    append_log "Filtering Results: $filtered_files_count / $total_iterations ($((filtered_files_count * 100 / total_iterations))%)"

    echo
}



# Check if the CSV file exists
if [ ! -e "$csv_file" ]; then
    # If the CSV file doesn't exist, create it and add the header
    echo "$header" > "$csv_file"
fi

# Call the functions
perform_terminus_scoring
perform_sesd_calculation
perform_substructure_counting
# filter_files
