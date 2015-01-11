
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
    - blastdb_exe - path to eht makeblastdb executable
    """
    cmdlines = []
    for idx, fname in enumerate(filenames[:-1]):
        cmdlines.extend([construct_makeblastdb_cmdline(fname, blastdb_exe) for 
                         fname in filenames[idx+1:]])
    return cmdlines


# Generate single makeblastdb command line
def construct_makeblastdb_cmdline(filename,
                                  blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastdb_exe - path to eht makeblastdb executable
    """
    title = os.path.splitext(os.path.split(filename)[-1])[0] + '-fragments'
    return "{0} -dbtype nucl -in {1} -title {2}".format(blastdb_exe,
                                                        filename, title)
