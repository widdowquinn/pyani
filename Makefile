# Makefile
#
# Runs tests for pyani module

DATA=tests/test_ani_data
OUT_B=tests/test_ANIb_output
OUT_BLASTALL=tests/test_ANIblastall_output
OUT_M=tests/test_ANIm_output
OUT_TETRA=tests/test_TETRA_output
CLASSES=$(DATA)/classes.tab
LABELS=$(DATA)/labels.tab

all : clean test

ANIm :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_M) \
	  -m ANIm --classes $(CLASSES) --labels $(LABELS) -g -v

ANIb :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_B) \
	  -m ANIb --classes $(CLASSES) --labels $(LABELS) -g -v

ANIblastall :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_BLASTALL) \
	  -m ANIblastall --classes $(CLASSES) --labels $(LABELS) -g -v

BLAST : ANIb ANIblastall

TETRA :
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_TETRA) \
	  -m TETRA --classes $(CLASSES) --labels $(LABELS) -g -v

clean :
	rm -rf $(OUT_M) $(OUT_B) $(OUT_BLASTALL) $(OUT_TETRA)

test : ANIm ANIb ANIblastall TETRA

