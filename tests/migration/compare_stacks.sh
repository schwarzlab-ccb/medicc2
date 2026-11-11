#!/usr/bin/env bash
# usage: tests/migration/compare_stacks.sh old|new    (run from the repo root, in the matching env)
#        tests/migration/compare_stacks.sh diff
set -u
OUT=${OUT:-/tmp/medicc_stack_compare}
EV="--events --chromosomes-bed default --regions-bed default --plot none --no-plot-tree"
runs=(
  "simple|examples/simple_example/simple_example.tsv $EV"
  "simple_j4|examples/simple_example/simple_example.tsv --n-cores 4 $EV"
  "testing|examples/testing_example/testing_example.tsv $EV"
  "testing_tcn|examples/testing_example/testing_example.tsv --total-copy-numbers --input-allele-columns cn_a $EV"
  "testing_nowgd|examples/testing_example/testing_example.tsv --no-wgd $EV"
  "testing_x2|examples/testing_example/testing_example.tsv --wgd-x2 $EV"
  "wgd|examples/WGD_example/WGD_example.tsv $EV"
  "ov|examples/OV03-04/OV03-04_descr.txt -i fasta --normal-name OV03-04_diploid $EV"
)
for p in PTX004 PTX005 PTX006 PTX007 PTX008 PTX009 PTX010 PTX011 PTX012 PTX013; do
  runs+=("g_$p|examples/gundem_et_al_2015/${p}_input_df.tsv $EV")
done
LEGACY=${LEGACY:-medicc2}          # becomes medicc2-legacy after Step 44
if [ "$1" = diff ]; then
  diff -r -x '*.png' -x '*.pdf' -x '*.log' $OUT/old $OUT/new && echo "PASS: identical outputs"; exit
fi
mkdir -p $OUT/$1
for r in "${runs[@]}"; do
  name=${r%%|*}; args=${r#*|}
  python $LEGACY $args $OUT/$1/$name > $OUT/$1/$name.log 2>&1 || echo "FAIL $name (see $OUT/$1/$name.log)"
done
echo "done: $OUT/$1"
