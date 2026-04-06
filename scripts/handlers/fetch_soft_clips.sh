#!/bin/bash
# Purpose: Fetch soft clip rasters (.tif files) from HPC to local directory
# Usage: ./fetch_soft_clips.sh
#
# Hardcoded paths:
#   HPC source: mccleary:/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/soft_clips
#   Species list: ~/Downloads/Kass1_species.txt
#   Local destination: ~/clim_risk_phylosdm/analysis/soft_clips

set -e  # Exit on error

# Hardcoded paths
HPC_PATH="mccleary:/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/soft_clips"
SPECIES_LIST="$HOME/Downloads/Kass1_species.txt"
LOCAL_DIR="$HOME/clim_risk_phylosdm/analysis/soft_clips"

# Validate species list file exists
if [[ ! -f "$SPECIES_LIST" ]]; then
    echo "Error: Species list file not found: $SPECIES_LIST"
    exit 1
fi

# Create local directory if it doesn't exist
mkdir -p "$LOCAL_DIR"

# Count total species
TOTAL=$(wc -l < "$SPECIES_LIST" | tr -d ' ')
echo "Fetching soft clip rasters for $TOTAL species from $HPC_PATH to $LOCAL_DIR"

# Read species list and fetch each tif
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
    TIF_FILE="${species}_soft_clip.tif"
    
    echo "[$COUNT/$TOTAL] Fetching $TIF_FILE..."
    
    if scp "${HPC_PATH}/${TIF_FILE}" "${LOCAL_DIR}/" 2>/dev/null; then
        echo "  ✓ Downloaded $TIF_FILE"
    else
        echo "  ✗ Failed to download $TIF_FILE"
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
echo "  Soft clip rasters saved to: $LOCAL_DIR"
echo "========================================"
