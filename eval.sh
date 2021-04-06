for file in tests/netlib/problem_files/*
do
   cargo flamegraph -o metrics/flamegraph/"$file".svg --bin rust-lp --features="binaries" -- "$file"
done