#!/bin/bash

# Define paths
BASE_DIR="/xdisk/twheeler/nsontakke/ATB_Analysis_0725"
COMMANDS_FILE="${BASE_DIR}/metabolic_commands.txt"
BATCH_DIR="${BASE_DIR}/metabolic_batches"

# Parameters
TOTAL_COMMANDS=21392
NUM_BATCHES=100

# Create batch directory
mkdir -p "$BATCH_DIR"

# Calculate commands per batch
COMMANDS_PER_BATCH=$((TOTAL_COMMANDS / NUM_BATCHES))  # 213
REMAINDER=$((TOTAL_COMMANDS % NUM_BATCHES))            # 92

echo "Creating batch command files..."
echo "Total commands: $TOTAL_COMMANDS"
echo "Number of batches: $NUM_BATCHES"
echo "Commands per batch: ~$COMMANDS_PER_BATCH"
echo "Output directory: $BATCH_DIR"
echo ""

# Create individual batch files
for BATCH_ID in $(seq 1 $NUM_BATCHES); do
    # Calculate start and end line for this batch
    if [ $BATCH_ID -le $REMAINDER ]; then
        # First batches get one extra command to distribute the remainder
        COMMANDS_THIS_BATCH=$((COMMANDS_PER_BATCH + 1))
        START_LINE=$(( (BATCH_ID - 1) * COMMANDS_THIS_BATCH + 1 ))
    else
        COMMANDS_THIS_BATCH=$COMMANDS_PER_BATCH
        START_LINE=$(( REMAINDER * (COMMANDS_PER_BATCH + 1) + (BATCH_ID - REMAINDER - 1) * COMMANDS_PER_BATCH + 1 ))
    fi
    
    END_LINE=$((START_LINE + COMMANDS_THIS_BATCH - 1))
    
    # Create batch file
    BATCH_FILE="${BATCH_DIR}/batch_${BATCH_ID}.txt"
    
    # Extract commands for this batch
    sed -n "${START_LINE},${END_LINE}p" "$COMMANDS_FILE" > "$BATCH_FILE"
    
    # Verify
    LINES_IN_FILE=$(wc -l < "$BATCH_FILE")
    
    echo "Batch $BATCH_ID: lines $START_LINE-$END_LINE ($LINES_IN_FILE commands) -> $BATCH_FILE"
done

echo ""
echo "Batch files created successfully!"
echo "Total batch files: $(ls -1 $BATCH_DIR/batch_*.txt | wc -l)"
echo ""
echo "Verify:"
echo "  Total lines in all batches: $(cat $BATCH_DIR/batch_*.txt | wc -l)"
