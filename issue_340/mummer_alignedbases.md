# How does `mummer` calculate `AlignedBases`?

`AlignedBases` is output from `dnadiff` and is (we think) meant to represent the number of non-redundant bases in the reference genome aligned to the query genome and vice versa. For example...

```bash
% dnadiff AF_bug_v2/MGV-GENOME-0357962.fna AF_bug_v2/MGV-GENOME-0358017.fna
% head -n 14 out.report
/Users/lpritc/Development/pyani/issue_340/donovan_test/AF_bug_v2/MGV-GENOME-0357962.fna /Users/lpritc/Development/pyani/issue_340/donovan_test/AF_bug_v2/MGV-GENOME-0358017.fna
NUCMER

                               [REF]                [QRY]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                     87285                87353
AlignedBases          87285(100.00%)       87353(100.00%)
UnalignedBases              0(0.00%)             0(0.00%)
```

Here, 87285 bases are aligned in the reference, 87353 in the query.

What does `mummer` do to calculate that value?

## `dnadiff.pl`

### Calling `mummer` programs to generate alignment output

In `dnadiff.pl`, some `mummer` programs are called (ll.124-129):

```perl
RunAlignment() unless defined($OPT_DeltaFile);
RunFilter();
RunCoords();
RunSNPs();
RunDiff();
MakeReport();
```

`RunAlignment()` uses `nucmer` with the command template:

```bash
$NUCMER --maxmatch -p $OPT_Prefix $OPT_RefFile $OPT_QryFile
```

This is a standard `nucmer` run, using the `--maxmatch` option, which uses all anchor matches regardless of their uniqueness.

`RunFilter` uses `delta-filter` with command template:

```bash
$DELTA_FILTER -1 $OPT_DeltaFile > $OPT_DeltaFile1
$DELTA_FILTER -m $OPT_DeltaFile > $OPT_DeltaFileM
```

This generates two delta files, one containing only 1-to-1 alignments (`.1delta`), and the other contains many-to-many alignments (`.mdelta`).

The script then identifies SNPs, but this is not strictly relevant to our analysis:

```bash
$SHOW_SNPS -rlTHC $OPT_DeltaFile1 > $OPT_SnpsFile
```

and calculates a difference file for each of the reference and query sequences (`.rdiff`, `.qdiff`):

```bash
$SHOW_DIFF -rH $OPT_DeltaFileM > $OPT_DiffRFile
$SHOW_DIFF -qH $OPT_DeltaFileM > $OPT_DiffQFile
```

The options suppress header output and specify which sequence is reported on.

The remainder of the script then processes these outputs to generate the report.

### Processing alignment output

Relevant variables in the script are (considering only the reference sequence data, at first):

- `rnSeqs`: number of sequences in the genome
- `rnBases`: number of bases in the genome
- `rnABases`: number of aligned bases in the genome; this is what gets reported as `AlignedBases` and is the focus of our analysis of the code.

The name of each sequence in the genome, and the corresponding length, is read into a hash table (l.252) analogous to a Python dictionary, keyed by sequence ID:

```perl
FastaSizes($OPT_RefFile, \%refs);
```

The total number of bases in the genome is derived from this hash table:

```perl
foreach ( values(%refs) ) {
      $rnSeqs++;
      $rnBases += $_;
  }
```

which increments the count of sequences, as well as the sum of sequence lengths (total bases).

The script opens the many-to-many coordinates file (`.mcoords`), which has the form:

```text
1	24024	63330	87353	24024	24024	99.98	87285	87353	27.52	27.50	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
23884	24176	1	293	293	293	100.00	87285	87353	0.34	0.34	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24107	87285	176	63368	63179	63193	99.92	87285	87353	72.38	72.34	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
```

This has 13 columns. Each line represents a single alignment. `dnadiff.pl` processes only some of these columns.

The number of alignments `rqnAlignsM++` is incremented with each line (total number of alignments, l.274). The sum of alignment region lengths is incremented with column 4 (`$rSumLenM += $A[4];`, l.275). The sums of match identities (`$rqSumIdyM += ($A[6] / 100.0) * ($A[4] + $A[5]);`) and match lengths (`$rqSumLenM += ($A[4] + $A[5]);`) are gathered.

**NOTE:** The alignments where columns 4 and 5 differ in value are those where we expect our problems to arise as this indicates at least one indel (though be aware that a deletion in the reference may be offset by a deletion in the query, also - this might in some cases make it look like the match is exact on both sequences when it is not).

Next, for this row/alignment in the `.mcoords` file, `dnadiff.pl` checks to see if this sequence from the genome has been seen before. If it _has **not**_, then it increments the number of aligned sequences `rnASeqs`, and updates the number of aligned bases `rnABases` (ll.280-285):

```perl
#-- If new ID, add to sequence and base count
if ( $refs{$A[11]} > 0 ) {
    $rnASeqs++;
    $rnABases += $refs{$A[11]};
    $refs{$A[11]} *= -1; # If ref has alignment, length will be -neg
}
```

It also updates breakpoint counts, but this is not strictly relevant to us.

Next, the difference files for the genome is opened (`.rdiff`). This has format:

```text
MGV_MGV-GENOME-0357962	JMP	24025	23883	-141
MGV_MGV-GENOME-0357962	GAP	24177	24106	-70	-118	48
```

The number of gaps is in column 4, and this is held in `$gap`; initially `$ins` (insertions) is set to this value also (l.332):

```perl
my $gap = $A[4];
my $ins = $gap;
```

If the `.rdiff` line is a `GAP` (column 1), the number of insertions (column 6) is assigned to `$ins`, **if that number is more positive than the value in `$gap`**. Then, if columns 4 and 5 have negative values, but column 6 is positive, the number of insertions (`$rnTIns`) is incremented, and the total length of insertions (`$rSumTIns`) is extended by the value in column 6.

If the `.rdiff` line is not a `DUP` then the number of aligned bases `$rnABases` has the value of `$gap` subtracted from it.

If `$ins` is positive (if there is an insertion), the number of insertions `$rnIns` is incremented, and the total insertion size (`$rSumIns`) is updated.

```perl
#-- Add to tandem insertion counts
if ( $A[1] eq "GAP" ) {
    scalar(@A) == 7
        or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
    $ins = $A[6] if ( $A[6] > $gap );
    if ( $A[4] <= 0 && $A[5] <= 0 && $A[6] > 0 ) {
        $rnTIns++;
        $rSumTIns += $A[6];
     }
}

#-- Remove unaligned sequence from count
if ( $A[1] ne "DUP" ) {
  $rnABases -= $gap if ( $gap > 0 );
}

#-- Add to insertion count
if ( $ins > 0 ) {
    $rnIns++;
    $rSumIns += $ins;
}
```

The script then deals with SNPs, which aren't important here, and goes on to calculate the number of aligned bases in ll.460-465:

```perl
    printf $fho "%-15s %20s %20s\n",
    "AlignedBases",
    ( sprintf "%10d(%.4f%%)",
      $rnABases, ($rnBases ? $rnABases / $rnBases * 100.0 : 0) ),
    ( sprintf "%10d(%.4f%%)",
      $qnABases, ($qnBases ? $qnABases / $qnBases * 100.0 : 0) );
```

to generate the output:

```text
AlignedBases          87285(100.00%)       87353(100.00%)
```

Here, `87285(100.00%)` is the reference genome aligned bases count. The `87285` is the value held in `$rnABases`. If the genome has any bases (`$rnBases`) it calculates the fraction of those which are aligned bases (`? $rnABases / $rnBases * 100.0`) else it reports zero (`: 0`).

By this point, the value in `$rnABases` represents:

- the sum of total lengths of all aligned (reference) sequences (ll.280-285)
- minus the number of gaps indicated for that sequence in the `.rdiff` file

## How is the `.rdiff` file generated?

The `.rdiff` file is generated by the `show-diff` `mummer` program, which has [documentation](https://mummer4.github.io/manual/manual.html):

> This program classifies alignment breakpoints for the quantification of macroscopic differences between two genomes. It takes a standard, unfiltered delta file as input, determines the best mapping between the two sequence sets, and reports on the breaks in that mapping.

```text
OUTPUT:
stdout  Classified breakpoints are output one per line with
        the following types and column definitions. The first
        five columns of every row are seq ID, feature type,
        feature start, feature end, and feature length.

Feature Columns

IDR GAP gap-start gap-end gap-length-R gap-length-Q gap-diff
IDR DUP dup-start dup-end dup-length
IDR BRK gap-start gap-end gap-length
IDR JMP gap-start gap-end gap-length
IDR INV gap-start gap-end gap-length
IDR SEQ gap-start gap-end gap-length prev-sequence next-sequence
```

## What is happening in our example?

### Raw alignment (`.delta`)

The example we're working with is an alignment of two viruses, `MGV-GENOME-0357962` and `MGV-GENOME-0358017`. Running `dnadiff.pl` generates an intial `.delta` file (`out.delta`). This is the primary output and all other actions use this data to generate human-readable/other output.

```text
/Users/lpritc/Development/pyani/issue_340/donovan_test/AF_bug_v2/MGV-GENOME-0357962.fna /Users/lpritc/Development/pyani/issue_340/donovan_test/AF_bug_v2/MGV-GENOME-0358017.fna
NUCMER
>MGV_MGV-GENOME-0357962 MGV_MGV-GENOME-0358017 87285 87353
1 24024 63330 87353 5 5 0
0
23884 24176 1 293 0 0 0
0
24107 87285 176 63368 51 51 0
-121
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-10416
-1
0
```

As per the documentation, this file gives us a lot of information.

- The header line indicates the input files (reference then query)
- The next lines indicates whether `NUCMER` or `PROMER` was used
- Then there is an alignment header, which states the aligned sequences, and the lengths of each:

```text
>MGV_MGV-GENOME-0357962 MGV_MGV-GENOME-0358017 87285 87353
```

There then follows a series of alignment data lines. Each new alignment begins with a header. Here, the first is:

```text
1 24024 63330 87353 5 5 0
```

which tells us that the alignment runs from base `1` to base `24024` in the reference sequences, and from base `63330` to base `87353` in the query sequence. The next three numbers are the number of "errors" (non-identities + indels, here `5`), the nunber of  "similarity errors" (non-positive match scores, here `5`), and stop codons (this is a DNA alignment, so this is `0`).

Each aligned region description ends with `0`, and for the first two aligned regions, there is no other information:

```text
1 24024 63330 87353 5 5 0
0
23884 24176 1 293 0 0 0
0
```

The third and final aligned region however has indels. The header indicates an alignment with `51`` non-identities and indels, and then data describing the match:

```text
24107 87285 176 63368 51 51 0
-121
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-10416
-1
0
```

The numbers following the header indicate the distance to the next indel in the alignment. A positive number means the indel is an insertion in the reference; a negative number means the indel is a deletion in the reference. Here we interpret this as:

- 121 aligned bases, then 12 deletions in the reference, then 10416 aligned bases, followed by two deletions in the reference.
  - This give us a total of 14 deletions in the reference alignment.

The overall picture this presents is the following three alignments:

```text
R      1 =========== 24024 (24024 bases)
Q  63330 =========== 87353 (24024 bases)

R  23884 =========== 24176 (293 bases)
Q      1 =========== 293   (293 bases)

R  24107 =========== 87285 (63179 bases, but 14 deletions)
Q    176 =========== 63368 (63193 bases)
```

It would be straightforward to calculate the reference genome's total number of aligned bases: `24024 + 293 + 63179 = 87496`, and this corresponds to the `TotalLength` value in `out.report` (as we would expect).

But these alignments overlap in the reference:

```text
1 ============ 24024
     23884 ============== 24176
                 24107 =========== 87285
```

So we need to subtract the duplicated overlapping stretches in each case: (`24024 - 23884 + 1 = 141`) and (`24176 - 24107 + 1 = 70`), for a total of `211` overlapping bases. The total number of aligned bases in the reference is then `87496 - 211 = 87275`, which is the expected number of `AlignedBases` in `out.report`.

But what of the query sequence? That doesn't have deletions, but it does have _insertions_ (relative to the reference). The total number of bases in the alignment is `24024 + 293 + 63193 = 87510`, as in `out.report`. The overlap sizes are (`293 - 176 + 1 = 118`) and (`63368 - 63330 + 1 = 39`) for a total of `157` bases, for total length of `87510 - 157 = 87353`, which matches `AlignedBases` in the `out.report`.

**But wait!** What about those 14 deletions in the reference genome? These correspond to bases in the query that have no counterpart in the reference genome. So are they aligned bases, if they do not align to the reference? I think they are not.

We have reproduced the `dnadiff.pl` `out.report` values, but we have not accounted for the insertions in the query/deletions in the reference! There are 14 of these. Why have they not been considered when calculating the number of aligned bases in the `out.report` file?

Checking the `out.snps` file, we can see the locations of the indels:

```text
23190	T	G	86519	835	835	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	296	0	296	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	G	297	0	297	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	298	0	298	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	299	0	299	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	G	300	0	300	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	301	0	301	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	G	302	0	302	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	G	303	0	303	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	304	0	304	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	305	0	305	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	G	306	0	306	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
24226	.	A	307	0	307	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
27841	C	T	3922	3615	3922	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
34641	.	A	10723	0	10723	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
34641	.	A	10724	0	10724	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
34782	C	T	10865	141	10865	87285	87353	1	1	MGV_MGV-GENOME-0357962	MGV_MGV-GENOME-0358017
```

The deletions are at bases 24226 and 34641 in the reference sequence, and insertions at bases 296 and 10723 in the query. Do they lie within one of the overlaps (and so can maybe be discounted)? A check against our alignments:

```text
R      1 =========== 24024 (24024 bases)
Q  63330 =========== 87353 (24024 bases)

R  23884 =========== 24176 (293 bases)
Q      1 =========== 293   (293 bases)

R  24107 =========== 87285 (63179 bases, but 14 deletions)
Q    176 =========== 63368 (63193 bases)
```

suggests that this is not the case. The deletions like in the long, non-overlapping section of the third alignment for both reference and query, and not in any of the other alignments.

A true count of aligned bases should take into account the failure to align bases due to indels. If we did so, we would adjust the count of aligned bases for the query to be `87353 - 14 = 87339`, for the 14 inserted bases in the query.

In practice, if we were to account for this, we would have to identify indels in overlapping regions. and handle them specifically to account for cases where the indel is not present in all of the overlap alignments. 

In our example here, the overlaps don't quite match up, either (`o[0-9]*` indicates overlap length):

```text
R      1 == 24024 (o141) 23884 == 24176 (o70) 24107 == 87285
Q  63330 == 87353 (o0)       1 == 293   (o118)  176 == 63368 (o39)
```

which means that a different number of bases are "duplicated" and align in the overlap for the reference and query. This may correspond to overlaps where there are more than two tandem repeats of similar short sequence, but different numbers of these repeats are identified in each alignment (or maybe to something else - it seems hard to word back from the numbers alone).

## What about many-to-many and one-to-one filtering?

