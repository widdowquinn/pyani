
import os
import pyani_config

from Bio import SeqIO

# Divide input FASTA sequences into fragments
def fragment_FASTA_files(infiles, outdirname, fragsize):
    """Chops sequences of the passed files into fragments, returns filenames.

    - infiles - paths to each input sequence file
    - outdirname - path to output directory
    - fragsize - the size of sequence fragments

    Takes every sequence from every file in infiles, and splits them into
    consecutive fragments of length fragsize, (with any trailing sequences
    being included, even if shorter than fragsize), and writes the resulting
    set of sequences to a file with the same name in the output directory.
    All fragments are named consecutively and uniquely (within a file) as
    fragNNNNN. Sequence description fields are retained.
    """
    outfnames = []
    for fname in infiles:
        outstem = os.path.splitext(os.path.split(fname)[-1])[0]
        outfname = os.path.join(outdirname, outstem) + '.fasta'
        outseqs = []
        idx, count = 0, 0
        for seq in SeqIO.parse(fname, 'fasta'):
            while idx < len(seq):
                count += 1
                newseq = seq[idx:idx+fragsize]
                newseq.id = "frag%05d" % count
                outseqs.append(newseq)
                idx += fragsize
        outfnames.append(outfname)
        SeqIO.write(outseqs, outfname, 'fasta')
    return outfnames


# Generate list of makeblastdb command lines from passed filenames
def generate_blastdb_commands(filenames, outdir,
                              blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT):
    """Return a list of makeblastdb command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - blastdb_exe - path to the makeblastdb executable
    """
    cmdlines = [construct_makeblastdb_cmdline(fname, blastdb_exe) for 
                fname in filenames]
    return cmdlines


# Generate single makeblastdb command line
def construct_makeblastdb_cmdline(filename,
                                  blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastdb_exe - path to the makeblastdb executable
    """
    title = os.path.splitext(os.path.split(filename)[-1])[0] + '-fragments'
    return "{0} -dbtype nucl -in {1} -title {2}".format(blastdb_exe,
                                                        filename, title)

# Generate list of BLASTN command lines from passed filenames
def generate_blastn_commands(filenames, outdir,
                             blastn_exe=pyani_config.BLASTN_DEFAULT):
    """Return a list of makeblastdb command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - blastn_exe - path to BLASTN executable
    """
    cmdlines = []
    for idx, fname1 in enumerate(filenames[:-1]):
        cmdlines.extend([construct_makeblastn_cmdline(fname1, fname2, outdir,
                                                      blastn_exe) for 
                         fname2 in filenames[idx+1:]])
    return cmdlines


# Generate single BLASTN command line
def construct_makeblastn_cmdline(fname1, fname2, outdir,
                                 blastn_exe=pyani_config.BLASTN_DEFAULT):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastn_exe - path to BLASTN executable
    """
    fstem1 = os.path.splitext(os.path.split(fname1)[-1])[0]
    fstem2 = os.path.splitext(os.path.split(fname2)[-1])[0]
    prefix = os.path.join(outdir, "%s_vs_%s" % (fstem1, fstem2))
    cmd = "{0} -out {1}.blast_tab -query {2} -db {3} " +\
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 " +\
        "-max_target_seqs 1 -outfmt '6 qseqid sseqid length mismatch " +\
        "pident nident qlen slen qstart qend sstart send positive " +\
        "ppos gaps' -task blastn"
    return cmd.format(blastn_exe, prefix, fname1, fname2) 


# Process pairwise BLASTN output
def process_blastn(outdirname, org_lengths):
    raise NotImplementedError
