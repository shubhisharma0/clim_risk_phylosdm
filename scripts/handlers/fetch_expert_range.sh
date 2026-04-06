#!/bin/bash
# Purpose: Fetch expert range maps (.gpkg files) from HPC to local directory
# Usage: ./fetch_expert_range.sh <hpc_path> <species_list_file> <local_dir>
#
# Arguments:
#   hpc_path          - Remote HPC path where .gpkg files are stored (e.g., mccleary:/path/to/ranges)
#   species_list_file - Path to text file containing species names (one per line) 
#   local_dir         - Local directory to store the downloaded .gpkg files
#
# Example:
#   ./fetch_expert_range.sh mccleary:/vast/palmer/pi/jetz/ss4224/expert_ranges species_list.txt ./expert_ranges

set -e  # Exit on error

# Check arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <hpc_path> <species_list_file> <local_dir>"
    echo "  hpc_path          - Remote HPC path (e.g., mccleary:/path/to/ranges)"
    echo "  species_list_file - Text file with species names (one per line)"
    echo "  local_dir         - Local directory to store .gpkg files"
    exit 1
fi

HPC_PATH="$1"
SPECIES_LIST="$2"
LOCAL_DIR="$3"

# Validate species list file exists
if [[ ! -f "$SPECIES_LIST" ]]; then
    echo "Error: Species list file not found: $SPECIES_LIST"
    exit 1
fi

# Create local directory if it doesn't exist
mkdir -p "$LOCAL_DIR"

# Count total species
TOTAL=$(wc -l < "$SPECIES_LIST" | tr -d ' ')
echo "Fetching expert range maps for $TOTAL species from $HPC_PATH to $LOCAL_DIR"

# Read species list and fetch each gpkg
COUNT=0
FAILED=0
LINENUM=0
while IFS= read -r species || [[ -n "$species" ]]; do
    LINENUM=$((LINENUM + 1))
    
    # Skip empty lines and comments
    [[ -z "$species" || "$species" =~ ^# ]] && continue
    
    # Strip quotes and trim whitespace
    species=$(echo "$species" | tr -d '"' | xargs)
    
    # Skip header row (if it contains "species" as column name)
    [[ "$species" == "species" || "$species" == "x" ]] && continue
    
    COUNT=$((COUNT + 1))
    GPKG_FILE="${species}.gpkg"
    
    echo "[$COUNT/$TOTAL] Fetching $GPKG_FILE..."
    
    if scp "${HPC_PATH}/${GPKG_FILE}" "${LOCAL_DIR}/" 2>/dev/null; then
        echo "  ✓ Downloaded $GPKG_FILE"
    else
        echo "  ✗ Failed to download $GPKG_FILE"
        FAILED=$((FAILED + 1))
    fi
done < "$SPECIES_LIST"

# Summary
echo ""
echo "========================================"
echo "Fetch complete!"
echo "  Successfully downloaded: $((COUNT - FAILED))/$COUNT"
if [[ $FAILED -gt 0 ]]; then
    echo "  Failed: $FAILED"
fi
echo "  Expert range maps saved to: $LOCAL_DIR"
echo "========================================"
