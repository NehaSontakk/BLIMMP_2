#!/usr/bin/env bash
# One-time setup: install BLIMMP Python dependencies into blimmp-work conda env
# and verify that all reference data files are in place.
#
# Run interactively on puma (not via sbatch):
#   bash setup_blimmp_env.sh
set -eo pipefail

source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate blimmp-work

echo "=== Installing BLIMMP Python dependencies into blimmp-work ==="
conda install -n blimmp-work -c conda-forge numba scipy pandas numpy -y

echo ""
echo "=== Verifying BLIMMP reference data ==="
BLIMMP_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/BLIMMP_8Sep2026/BLIMMP"

REQUIRED=(
    "ATB Frequency"
    "ONE_HOP_NEIGHBOR_DATA-2"
    "TWO_HOP_NEIGHBOR_DATA-2"
    "MODULE_ALL_NEIGHBOR_DATA"
    "KEGG_Module_Equations_18AUG26.json"
    "KEGG_Graphs_Generated_13AUG2026"
    "ko_list.txt"
    "module_freq.txt"
    "module_ko_reaction.json"
    "kegg_bacteria_modules.json"
    "ko_reaction.list"
    "CoOccurrence_AllVsAll_domain_level_priors.json"
    "BLIMMP_Version_26Sep2026.py"
)

MISSING=0
for rel in "${REQUIRED[@]}"; do
    full="$BLIMMP_DIR/$rel"
    if [[ -e "$full" ]]; then
        printf "  OK  : %s\n" "$rel"
    else
        printf "  MISS: %s\n" "$rel"
        MISSING=$(( MISSING + 1 ))
    fi
done

echo ""
if (( MISSING > 0 )); then
    echo "MISSING $MISSING item(s) from $BLIMMP_DIR"
    exit 1
else
    echo "All reference data found."
fi

echo ""
echo "=== Quick Python import test ==="
python - <<'PYEOF'
import sys
try:
    import numba; import scipy.stats; import pandas; import numpy
    print(f"  numba  {numba.__version__}")
    print(f"  scipy  {scipy.__version__}")
    print(f"  pandas {pandas.__version__}")
    print(f"  numpy  {numpy.__version__}")
    print("All imports OK.")
except ImportError as e:
    print(f"Import error: {e}"); sys.exit(1)
PYEOF

echo ""
echo "=== Test wrapper import ==="
WRAPPER="$(dirname "$0")/blimmp_puma_wrapper.py"
if [[ -f "$WRAPPER" ]]; then
    python - "$WRAPPER" <<'PYEOF'
import sys, importlib.util
spec = importlib.util.spec_from_file_location("w", sys.argv[1])
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("Wrapper import OK. BlimmpPipeline accessible:", hasattr(m.blimmp, 'BlimmpPipeline'))
PYEOF
else
    echo "blimmp_puma_wrapper.py not found next to this script — copy it to $(dirname "$0")/"
fi

echo ""
echo "Setup complete. You can now submit:"
echo "  sbatch run_blimmp_spliced.sh"