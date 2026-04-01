#!/bin/bash
#SBATCH --job-name=new_csv_bin_job
#SBATCH --output=logs/new_csv_bin_%j.out
#SBATCH --error=logs/new_csv_bin_%j.err
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --partition=all

# New CSV Bin SLURM Job Submitter
# Usage: ./slurm_new_csv_bin input_files
# 
# Where input_files folder contains:
#   - input_table.csv (bin definitions)
#   - materials.csv (material properties)
#   - mesh.py (mesh definition)
#   - scenario_*.py (scenario file - any name ending in .py except mesh.py and temperature_models.py)
#   - temperature_models.py (optional per-material temperature model overrides)

# Load modules and activate environment
module load IMAS
source /home/ITER/llealsa/miniconda3/etc/profile.d/conda.sh
conda activate PFC-TT
export PATH="/home/ITER/llealsa/miniconda3/envs/PFC-TT/bin:$PATH"

unset PYTHONPATH
export PYTHONNOUSERSITE=1

module unload SciPy-bundle        2>/dev/null
module unload Python-bundle-PyPI  2>/dev/null
module unload Python              2>/dev/null
module unload numpy               2>/dev/null
module unload mpi4py              2>/dev/null
module unload scifem              2>/dev/null

# Parse command line arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_files [bin_specification]"
    echo ""
    echo "Where input_files is a folder containing:"
    echo "  - input_table.csv"
    echo "  - materials.csv"
    echo "  - mesh.py"
    echo "  - scenario_*.py (or any .py file except mesh.py)"
    echo ""
    echo "Examples:"
    echo "  $0 input_files                  # Run all bins"
    echo "  $0 input_files \"0-4\"            # Run bins 0 to 4"
    echo "  $0 input_files \"0-4, 10-15\"     # Run bins 0-4 and 10-15"
    exit 1
fi

INPUT_DIR="$1"
BIN_SPEC="${@:2}"  # Everything after input_dir

# ---------- Resolve input directory ----------
# Accepts an absolute path, a relative path, or a bare folder name.
# For a bare folder name (no slashes), search both inside the current
# directory (PFC-Tritium-Transport) and one level above it.
# Errors out if the folder is found in both locations (ambiguous).
case "$INPUT_DIR" in
    /*)
        # Absolute path - use as-is
        if [ ! -d "$INPUT_DIR" ]; then
            echo "Error: Input directory '$INPUT_DIR' not found!"
            exit 1
        fi
        ;;
    */*)
        # Relative path with slashes - use as-is, then make absolute
        if [ ! -d "$INPUT_DIR" ]; then
            echo "Error: Input directory '$INPUT_DIR' not found!"
            exit 1
        fi
        INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
        ;;
    *)
        # Bare folder name - search in CWD and one level up
        _in_cwd=false
        _in_parent=false
        _cwd_abs=""
        _parent_abs=""

        if [ -d "./$INPUT_DIR" ]; then
            _in_cwd=true
            _cwd_abs="$(cd "./$INPUT_DIR" && pwd)"
        fi
        if [ -d "../$INPUT_DIR" ]; then
            _in_parent=true
            _parent_abs="$(cd "../$INPUT_DIR" && pwd)"
        fi

        # Avoid false positive when both resolve to the same directory
        if $_in_cwd && $_in_parent && [ "$_cwd_abs" != "$_parent_abs" ]; then
            echo "Error: Input directory '$INPUT_DIR' found in BOTH locations:"
            echo "  Inside PFC-TT:  $_cwd_abs"
            echo "  Outside PFC-TT: $_parent_abs"
            echo "  Provide a full or relative path to disambiguate."
            exit 1
        fi

        if $_in_cwd; then
            INPUT_DIR="$_cwd_abs"
        elif $_in_parent; then
            echo "  Resolved input directory (outside PFC-TT): $_parent_abs"
            INPUT_DIR="$_parent_abs"
        else
            echo "Error: Input directory '$INPUT_DIR' not found!"
            echo "  Searched: $(pwd)/$INPUT_DIR"
            echo "  Searched: $(cd .. && pwd)/$INPUT_DIR"
            exit 1
        fi
        ;;
esac

# Ensure INPUT_DIR is absolute (for paths that were already absolute above)
case "$INPUT_DIR" in
    /*) ;; # already absolute
    *)  INPUT_DIR="$(cd "$INPUT_DIR" && pwd)" ;;
esac

# Check for required files
if [ ! -f "$INPUT_DIR/input_table.csv" ]; then
    echo "Error: $INPUT_DIR/input_table.csv not found!"
    exit 1
fi

if [ ! -f "$INPUT_DIR/materials.csv" ]; then
    echo "Error: $INPUT_DIR/materials.csv not found!"
    exit 1
fi

if [ ! -f "$INPUT_DIR/mesh.py" ]; then
    echo "Error: $INPUT_DIR/mesh.py not found!"
    exit 1
fi

# Find the scenario file (any .py file except mesh.py and optional temperature_models.py)
SCENARIO_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.py" ! -name "mesh.py" ! -name "temperature_models.py" | head -1)
if [ -z "$SCENARIO_FILE" ] || [ ! -f "$SCENARIO_FILE" ]; then
    echo "Error: No scenario .py file found in '$INPUT_DIR' (excluding mesh.py and temperature_models.py)!"
    exit 1
fi

# Extract scenario name from filename (without .py extension)
SCENARIO_NAME=$(basename "$SCENARIO_FILE" .py)

CSV_FILE="$INPUT_DIR/input_table.csv"
SCENARIO_FOLDER="$INPUT_DIR"

echo "Input Configuration:"
echo "  Input folder: $INPUT_DIR"
echo "  CSV file: $CSV_FILE"
echo "  Materials file: $INPUT_DIR/materials.csv"
echo "  Mesh file: $INPUT_DIR/mesh.py"
echo "  Scenario file: $SCENARIO_FILE"
echo "  Scenario name: $SCENARIO_NAME"

# Function to expand bin specifications into individual bin IDs
expand_bin_spec() {
    local spec="$1"
    local bins=()
    
    # Handle comma-separated ranges and individual numbers
    spec=$(echo "$spec" | tr ',' ' ')
    
    for token in $spec; do
        token=$(echo "$token" | xargs)
        
        if [[ $token =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start=${BASH_REMATCH[1]}
            end=${BASH_REMATCH[2]}
            for ((i=start; i<=end; i++)); do
                bins+=($i)
            done
        elif [[ $token =~ ^[0-9]+$ ]]; then
            bins+=($token)
        else
            echo "Error: Invalid bin specification '$token'. Use format like '0-5', '10', or '0-5, 10, 15-20'"
            exit 1
        fi
    done
    
    printf '%s\n' "${bins[@]}" | sort -n | uniq
}

# Determine which sim IDs to run
# If a "Sim. ID" (or "sim_id" etc.) column exists in the CSV, use those
# values as the sim IDs.  Otherwise fall back to 1-based row numbers.
#
# When a BIN_SPEC is provided on the command line, those numbers refer to
# sim IDs (not row numbers).

# Helper: extract sim IDs from CSV header + data
get_sim_ids_from_csv() {
    local csv="$1"
    # Read the header line and normalise (lowercase, strip spaces/underscores/dots)
    local header
    header=$(head -1 "$csv")
    # Find the 1-based column index of a "sim id" column
    local col_idx
    col_idx=$(echo "$header" | awk -F',' '{
        for (i=1; i<=NF; i++) {
            col = tolower($i);
            gsub(/[ _.\-()]/, "", col);
            if (col == "simid" || col == "simulationid") { print i; exit }
        }
    }')
    if [ -n "$col_idx" ]; then
        # Extract the sim_id column values (skip header)
        tail -n +2 "$csv" | awk -F',' -v c="$col_idx" '{
            val = $c; gsub(/[ \t]/, "", val);
            if (val != "") print val
        }'
    else
        # No sim_id column → generate 1..NUM_ROWS
        local num_rows
        num_rows=$(awk 'END{print NR-1}' "$csv")
        seq 1 "$num_rows"
    fi
}

if [ -z "$BIN_SPEC" ]; then
    # Run ALL sim IDs from the CSV
    mapfile -t SIM_IDS_ARRAY < <(get_sim_ids_from_csv "$CSV_FILE")
    echo "  Sim ID source: CSV column (or 1-based rows)"
    echo "  Total sims: ${#SIM_IDS_ARRAY[@]}"
else
    BIN_IDS_OUTPUT=$(expand_bin_spec "$BIN_SPEC")
    SIM_IDS_ARRAY=($BIN_IDS_OUTPUT)
    echo "  Sim ID specification: $BIN_SPEC"
    echo "  Total sims: ${#SIM_IDS_ARRAY[@]}"
fi

# Create logs directory inside the input folder
LOGS_DIR="$INPUT_DIR/logs"
mkdir -p "$LOGS_DIR"

echo ""
echo "Submitting jobs to SLURM cluster..."
echo "=========================================="

# Loop over specified sim IDs
for sim_id in "${SIM_IDS_ARRAY[@]}"; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=sim_${sim_id}
#SBATCH --output=${LOGS_DIR}/sim_${sim_id}_%j.out
#SBATCH --error=${LOGS_DIR}/sim_${sim_id}_%j.err
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --partition=all

module load IMAS
source /home/ITER/llealsa/miniconda3/etc/profile.d/conda.sh
conda activate PFC-TT
export PATH="/home/ITER/llealsa/miniconda3/envs/PFC-TT/bin:\$PATH"

unset PYTHONPATH
export PYTHONNOUSERSITE=1

module unload SciPy-bundle        2>/dev/null
module unload Python-bundle-PyPI  2>/dev/null
module unload Python              2>/dev/null
module unload numpy               2>/dev/null
module unload mpi4py              2>/dev/null
module unload scifem              2>/dev/null

PYTHONNOUSERSITE=1 python -s run_on_cluster/run_new_csv_bin.py $sim_id $SCENARIO_FOLDER $SCENARIO_NAME $CSV_FILE --input-dir $INPUT_DIR

EOF
    echo "Submitted job for Sim ID: $sim_id"
done

echo ""
echo "=========================================="
echo "All jobs submitted!"
echo "  Input folder: $INPUT_DIR"
echo "  Scenario: $SCENARIO_NAME"
echo "  Jobs submitted: ${#SIM_IDS_ARRAY[@]}"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Cancel: scancel -u \$USER"
