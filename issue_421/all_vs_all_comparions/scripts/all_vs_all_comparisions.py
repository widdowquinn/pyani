from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import os
import time


def construct_nucmer_cmdline(
    fname1: Path,
    fname2: Path,
    outdir: Path = Path("."),
    nucmer_exe = Path("nucmer"),
    filter_exe = Path("delta-filter"),
    maxmatch: bool = False,
) -> Tuple[str, str]:
    """Return a tuple of corresponding NUCmer and delta-filter commands.

    :param fname1:  path to query FASTA file
    :param fname2:  path to subject FASTA file
    :param outdir:  path to output directory
    :param nucmer_exe:
    :param filter_exe:
    :param maxmatch:  Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead

    The split into a tuple was made necessary by changes to SGE/OGE.
    The delta-filter command must now be run as a dependency of the NUCmer
    command, and be wrapped in a Python script to capture STDOUT.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".
    """
    # Cast path strings to pathlib.Path for safety
    fname1, fname2 = [Path(fname1), Path(fname2)]

    # Compile commands
    # Nested output folders to avoid N^2 scaling in files-per-folder
    # Create folders incrementally (want an error if outdir does not exist)
    outsubdir = outdir / "nucmer_output"
    outsubdir.mkdir(exist_ok=True)
    outsubdir = outdir / "nucmer_output" / fname1.stem
    outsubdir.mkdir(exist_ok=True)
    outprefix = outsubdir / f"{fname1.stem}_vs_{fname2.stem}"
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = "{0} {1} -p {2} {3} {4}".format(
        nucmer_exe, mode, outprefix, fname1, fname2
    )
    # There's a subtle pathlib.Path issue, here. We must use string concatenation to add suffixes
    # to the outprefix files, as using path.with_suffix() instead can replace part of the filestem
    # in those cases where there is a period in the stem (this occurs frequently as it is part
    # of the NCBI notation for genome assembly versions)
    filtercmd = (
        f"delta_filter_wrapper.py {filter_exe} -1 {str(outprefix) + '.delta'} "
        f"{str(outprefix) + '.filter'}"
    )
    return (nucmercmd, filtercmd)




def generate_nucmer_commands(
    filenames: List[Path],
    outdir: Path = Path("."),
    nucmer_exe = Path("nucmer"),
    filter_exe = Path("delta-filter"),
    maxmatch: bool = False,
) -> Tuple[List, List]:
    """Return list of NUCmer command-lines for ANIm.

    :param filenames:  a list of paths to input FASTA files
    :param outdir:  path to output directory
    :param nucmer_exe:  location of the nucmer binary
    :param maxmatch:  Boolean flag indicating to use NUCmer's -maxmatch option

    The first element returned is a list of NUCmer commands, and the
    second a corresponding list of delta_filter_wrapper.py commands.
    The NUCmer commands should each be run before the corresponding
    delta-filter command.

    TODO: This return value needs to be reworked as a collection.

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    nucmer_cmdlines, delta_filter_cmdlines = [], []
    filenames = sorted(filenames)  # enforce ordering of filenames
    for idx, fname1 in enumerate(filenames[:-1]):
        for fname2 in filenames[idx + 1 :]:
            ncmd, dcmd = construct_nucmer_cmdline(
                fname1, fname2, outdir, nucmer_exe, filter_exe, maxmatch
            )
            nucmer_cmdlines.append(ncmd)
            delta_filter_cmdlines.append(dcmd)


            ncmd_rvs, dcmd_rvs = construct_nucmer_cmdline(
                fname2, fname1, outdir, nucmer_exe, filter_exe, maxmatch
            )
            nucmer_cmdlines.append(ncmd_rvs)
            delta_filter_cmdlines.append(dcmd_rvs)


    return (nucmer_cmdlines, delta_filter_cmdlines)





class Job(object):

    """Individual job to be run, with list of dependencies."""

    def __init__(self, name: str, command: str, queue: Optional[str] = None) -> None:
        """Instantiate a Job object.

        :param name:           String describing the job (uniquely)
        :param command:        String, the valid shell command to run the job
        :param queue:          String, the SGE queue under which the job shall run
        """
        self.name = name  # Unique name for the job
        self.queue = queue  # The SGE queue to run the job under
        self.command = command  # Command line to run for this job
        self.script = command
        self.scriptPath = None  # type: Optional[Any]
        self.dependencies = []  # type: List[Any]
        self.submitted = False  # type: bool
        self.finished = False  # type: int

    def add_dependency(self, job) -> None:
        """Add passed job to the dependency list for this Job.

        :param job:  Job to be added to the Job's dependency list

        This Job should not execute until all dependent jobs are completed.
        """
        self.dependencies.append(job)

    def remove_dependency(self, job) -> None:
        """Remove passed job from this Job's dependency list.

        :param job:     Job to be removed from the Job's dependency list
        """
        self.dependencies.remove(job)

    def wait(self, interval: float = 0.01) -> None:
        """Wait until the job finishes, and poll SGE on its status.

        :param interval:  float, number of seconds to wait before polling SGE
        """
        while not self.finished:
            time.sleep(interval)
            interval = min(2.0 * interval, 60)
            self.finished = os.system(f"qstat -j {self.name} > /dev/null")

# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(
    filenames: List[Path],
    outdir: Path = Path("."),
    nucmer_exe = Path("nucmer"),
    filter_exe = Path("delta-filter"),
    maxmatch: bool = False,
    jobprefix: str = "ANINUCmer",
):
    """Return list of Jobs describing NUCmer command-lines for ANIm.

    :param filenames:  Iterable, Paths to input FASTA files
    :param outdir:  str, path to output directory
    :param nucmer_exe:  str, location of the nucmer binary
    :param filter_exe:
    :param maxmatch:  Boolean flag indicating to use NUCmer's -maxmatch option
    :param jobprefix:

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    ncmds, fcmds = generate_nucmer_commands(
        filenames, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ncmd in enumerate(ncmds):
        print(ncmd)
        njob = Job(f"{jobprefix}_{idx:06d}-n", ncmd)
        fjob = Job(f"{jobprefix}_{idx:06d}-f", fcmds[idx])
        fjob.add_dependency(njob)
        joblist.append(fjob)


    return joblist



generate_nucmer_jobs(["../../symmetry_data/input/test_1/MGV-GENOME-0264574.fna", "../../symmetry_data/input/test_1/MGV-GENOME-0266457.fna"])