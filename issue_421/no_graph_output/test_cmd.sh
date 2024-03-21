# Removing current pyANI database, and creating new one
rm -rf .pyani
pyani createdb

# Running analysis on 2 streptomyces genomes
pyani fastani -i  streptomyces/input -o  streptomyces/output -l test_1.log \
   --name "test_1" --labels streptomyces/input/labels.txt --classes  streptomyces/input/classes.txt \
   --debug
# # Get reports
pyani report --runs -o  streptomyces/output --formats=stdout --run_results 1
# # Get matrices
pyani report -o  streptomyces/output --formats=stdout --run_matrices 1 --debug
pyani plot -o streptomyces/output --run_id 1 --formats svg --debug

#Running analysis on 2 viral genomes (donovan)
pyani anim -i donovan_test/input  -o  donovan_test/output -l test_2.log \
   --name "test_2" --labels donovan_test/input/labels.txt --classes  donovan_test/input/classes.txt \
   --debug
# Get reports
pyani report --runs -o  donovan_test/output --formats=stdout --run_results 2
# # Get matrices
pyani report -o  donovan_test/output --formats=stdout --run_matrices 2
pyani plot -o donovan_test/output --run_id 2 -v --formats pdf

#Running analysis on 2 viral genomes (donovan 2)
pyani anim -i donovan_test_2/input  -o  donovan_test_2/output -l test_3.log \
   --name "test_3" --labels donovan_test_2/input/labels.txt --classes  donovan_test_2/input/classes.txt \
   --debug
# Get reports
pyani report --runs -o  donovan_test_2/output --formats=stdout --run_results 3
# # Get matrices
pyani report -o  donovan_test_2/output --formats=stdout --run_matrices 3 --debug
pyani plot -o donovan_test_2/output --run_id 3 -v --formats pdf

