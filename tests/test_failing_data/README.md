# README.md - test_failing_data

This set of sequences should be a 'failing' data set, where at least one pair of genome sequences has no appreciable similarity to each other, generating zero values from alignment output that can be problematic in processing.

## Genomes

* `NC_000918.fna`: Aquifex aeolicus VF5 chromosome
* `NC_010473.fna`: Escherichia coli str. K-12 substr. DH10B chromosome
* `NC_013353.fna`: Escherichia coli O103:H2 str. 12009
* `NC_015216.fna`: Methanobacterium sp. AL-21 chromosome
* `NC_023044.fna`: Methanobacterium sp. MB1 complete sequence

## Command-line

```graph    
bin/average_nucleotide_identity.py -i tests/test_failing_data/ -o ./test_out -m ANIm -g --gformat png
```

## Expected error

```
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_015216_vs_NC_000918.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_010473_vs_NC_023044.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_000918_vs_NC_023044.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_010473_vs_NC_015216.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_013353_vs_NC_015216.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_010473_vs_NC_000918.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_013353_vs_NC_000918.delta is zero!
WARNING: Total alignment length reported in ./test_out/nucmer_output/NC_013353_vs_NC_023044.delta is zero!
ERROR: This is possibly due to a NUCmer comparison being too distant for use. Please consider using the --maxmatch option.
ERROR: This is alternatively due to NUCmer run failure, analysis will continue, but please investigate.
```