if [ ! -d metrics ]; then
  mkdir metrics;
  mkdir metrics/flamegraph;
fi

echo "STARTING NEW MEASUREMENT" >> metrics/timings.txt

for file in tests/netlib/problem_files/*
do
   start=$SECONDS
   name=$(basename $file)
   cargo flamegraph -o metrics/flamegraph/"$name".svg --bin rust-lp --features="binaries" -- "$file"
   duration=$((SECONDS - start))
   echo "$file completed in $duration seconds" >> metrics/timings.txt
done
echo " " >> metrics/timings.txt