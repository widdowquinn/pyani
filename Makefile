# Makefile
#
# Runs tests for pyani module

DATA=tests/test_ani_data
OUT_B=tests/test_ANIb_output
OUT_BLASTALL=tests/test_ANIblastall_output
OUT_M=tests/test_ANIm_output
OUT_TETRA=tests/test_TETRA_output

all : clean test

ANIb :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_M) -m ANIm -g -v

ANIm :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_B) -m ANIb -g -v 

ANIblastall :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_BLASTALL) \
	  -m ANIblastall -g -v

TETRA :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_TETRA) \
	  -m TETRA -g -v

clean :
	rm -rf $(OUT_M) $(OUT_B) $(OUT_BLASTALL) $(OUT_TETRA)

test : ANIm ANIb ANIblastall TETRA

