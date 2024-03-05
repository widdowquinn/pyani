rm -r .pyani
pyani createdb

# Run on real genomes
pyani anim -i symmetry_data/input/test_1 -o symmetry_data/output/test_1 -l test_1.log \
   --name "test_1" --labels symmetry_data/input/test_1/labels.txt --classes symmetry_data/input/test_1/classes.txt \
   --debug
mkdir -p 2024-03-05_test_real/
pyani report --run_results 1 --run_matrices 1 -o 2024-03-05_test_real/ --debug

# Run on synthetic genomes where we know the answer
# Two identical 100bp "genomes" except that one genome has a 2bp deletion within the sequence
# This gives a single alignment that runs through the deletion
pyani anim -i aln_length_issue/input -o aln_length_issue/output -l test_2.log \
   --name "test_2" --labels aln_length_issue/input/labels.txt --classes aln_length_issue/input/classes.txt \
   --debug
mkdir -p 2024-03-05_test_synth/   
pyani report --run_results 2 --run_matrices 2 -o 2024-03-05_test_synth/ --debug
