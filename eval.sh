#!/bin/bash

function clean_up {
  echo " " >> metrics/timings.txt
  >Cargo.toml
  sed '!d' temp.toml >> Cargo.toml
  rm temp.toml
  rm perf.data
  rm perf.data.old
}


control_c() {
  clean_up
  pkill -P $$
  exit
}


trap control_c SIGINT


if [ ! -d metrics ]; then
  mkdir metrics;
  mkdir metrics/flamegraph;
fi

cp Cargo.toml temp.toml
>Cargo.toml
sed '/^\[profile.release\].*/a debug=true' temp.toml >> Cargo.toml

echo "STARTING NEW MEASUREMENT" >> metrics/timings.txt

read -e -p "Enter netlib test prefix: " prefix

for file in tests/netlib/problem_files/"$prefix"*
do
   start=$SECONDS
   name=$(basename $file)
   cargo flamegraph -o metrics/flamegraph/"$name".svg --bin rust-lp --features="binaries" -- "$file"
   duration=$((SECONDS - start))
   echo "$file completed in $duration seconds" >> metrics/timings.txt
done

clean_up
