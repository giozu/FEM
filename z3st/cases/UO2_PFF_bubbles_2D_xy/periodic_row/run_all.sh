#!/bin/bash
export PATH=/home/baptiste/miniconda3/envs/z3st/bin:$PATH
export LD_LIBRARY_PATH=/home/baptiste/miniconda3/envs/z3st/lib:$LD_LIBRARY_PATH
CASE=/home/baptiste/z3st/z3st/cases/UO2_PFF_bubbles_2D_xy/periodic_row
for s in 400 200 120 90 72 66 60; do
  echo "=== starting s_$s at $(date) ==="
  cd "$CASE/s_$s"
  /home/baptiste/miniconda3/envs/z3st/bin/python -m z3st > log_z3st.md 2>&1
  ec=$?
  nsucc=$(grep -c SUCCESS log_z3st.md)
  echo "=== s_$s done at $(date): exit=$ec success_count=$nsucc ===" 
done
echo "ALL_DONE"
