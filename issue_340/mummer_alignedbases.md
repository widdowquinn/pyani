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
- minus the number of gaps indicated for that sequence in the `.mcoords` file