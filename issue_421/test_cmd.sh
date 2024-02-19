rm -rf ./pyani
pyani createdb
pyani anim -i symmetry_data/input/test_1 -o symmetry_data/output/test_1 -l test_1.log \
   --name "test_1" --labels symmetry_data/input/test_1/labels.txt --classes symmetry_data/input/test_1/classes.txt \
   --debug
