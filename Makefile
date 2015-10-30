# Makefile
#
# Runs tests for pyani module
#
# =====
# LAZY RUNNING: Don't redo calls to BLAST and MUMmer
#
# If we want to run tests without clobbering old data/redoing calculations
# then use the argument -e LAZY=1, e.g.
#
# make ANIm -e LAZY=1
#
# NOTE: if you do this, don't use 'make all', as it will clean out the 
#       necessary intermediate files; use 'make test' instead to run everything
# =====

DATA=tests/test_ani_data
OUT_B=tests/test_ANIb_output
OUT_BLASTALL=tests/test_ANIblastall_output
OUT_M=tests/test_ANIm_output
OUT_TETRA=tests/test_TETRA_output
CLASSES=$(DATA)/classes.tab
LABELS=$(DATA)/labels.tab

# Decide whether we're skipping generation of intermediate files
ifeq ($(LAZY), 1)
	LAZYBLAST=--force --noclobber --skip_blast
	LAZYMUMMER=--force --noclobber --skip_nucmer
	LAZYTETRA=--force --noclobber
else
	LAZYBLAST=
	LAZYMUMMER=
	LAZYTETRA=
endif

all : clean test

sge : clean test_sge

ANIm :
	@echo "Testing ANIm analysis..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_M) \
	  -m ANIm --classes $(CLASSES) --labels $(LABELS) -g -v \
	  $(LAZYMUMMER)
	@echo "ANIm analysis test done"
	@echo

ANIm_sge :
	@echo "Testing ANIm analysis with SGE..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_M) \
	  -m ANIm --classes $(CLASSES) --labels $(LABELS) -g -v \
	  --scheduler SGE $(LAZYMUMMER)
	@echo "ANIm analysis test done"
	@echo

ANIb :
	@echo "Testing ANIb analysis..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_B) \
	  -m ANIb --classes $(CLASSES) --labels $(LABELS) -g -v \
	  $(LAZYBLAST)
	@echo "ANIb analysis test done"
	@echo

ANIb_sge :
	@echo "Testing ANIb analysis with SGE..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_B) \
	  -m ANIb --classes $(CLASSES) --labels $(LABELS) -g -v \
	  --scheduler SGE $(LAZYBLAST)
	@echo "ANIb analysis test done"
	@echo

ANIblastall :
	@echo "Testing ANIblastall analysis..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_BLASTALL) \
	  -m ANIblastall --classes $(CLASSES) --labels $(LABELS) -g -v \
	  $(LAZYBLAST)
	@echo "ANIblastall analysis test done"
	@echo

ANIblastall_sge :
	@echo "Testing ANIblastall analysis with SGE..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_BLASTALL) \
	  -m ANIblastall --classes $(CLASSES) --labels $(LABELS) -g -v \
	  --scheduler SGE $(LAZYBLAST)
	@echo "ANIblastall analysis test done"
	@echo

BLAST : ANIb ANIblastall

TETRA :
	@echo "Testing TETRA analysis..."
	./average_nucleotide_identity.py -i $(DATA) -o $(OUT_TETRA) \
	  -m TETRA --classes $(CLASSES) --labels $(LABELS) -g -v \
	  $(LAZYTETRA)
	@echo "TETRA analysis test done"
	@echo

clean :
	rm -rf $(OUT_M) $(OUT_B) $(OUT_BLASTALL) $(OUT_TETRA)

test : ANIm ANIb ANIblastall TETRA

test_sge : ANIm_sge ANIb_sge ANIblastall_sge
