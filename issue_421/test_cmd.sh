# Removing current pyANI database, and creating new one
rm -rf .pyani
pyani createdb

#Running analysis on 2 fake genomes
pyani anim -i  aln_length_issue/input -o  aln_length_issue/output -l test_1.log \
   --name "test_1" --labels aln_length_issue/input/labels.txt --classes  aln_length_issue/input/classes.txt \
   --debug
# Get reports
pyani report --runs -o  aln_length_issue/output --formats=stdout --run_results 1
# # Get matrices
pyani report -o  aln_length_issue/output --formats=stdout --run_matrices 1 --debug
pyani plot -o aln_length_issue/output --run_id 1 -v --formats pdf

#Running analysis on 2 viral genomes (donovan)
pyani anim -i   -o  donovan_test/output -l test_2.log \
   --name "test_2" --labels donovan_test/input/labels.txt --classes  donovan_test/input/classes.txt \
   --debug
# Get reports
pyani report --runs -o  donovan_test/output --formats=stdout --run_results 2
# # Get matrices
pyani report -o  donovan_test/output --formats=stdout --run_matrices 2 --debug
pyani plot -o donovan_test/output --run_id 2 -v --formats pdf
