# -*- coding: utf-8 -*-
"""Module providing useful functions for manipulating pyani's SQLite3 db.

(c) The James Hutton Institute 2016-2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2016-2018 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np
import pandas as pd
import sqlite3

# SQL SCRIPTS
# ============
# The following is SQL for various database admin tasks, defined here to
# be out of the way when reading code.

# Create database tables
#
# genomes            - a row of data, one per genome
# runs               - a row of data, one per pyani run
# comparisons        - a row of data, one per pairwise comparison
# run_genomes        - table providing many:many comparisons
#
# The intention is that a run applies to some/all genomes from the genomes
# table, and that all the relevant pairwise comparisons/results are stored
# in the comparisons table.
#
# Information about the run (when it was run, what command/method, etc.) are
# stored in the runs table.
#
# All genomes (whether used or not) are described in the genomes table. An MD5
# hash is used to uniquely identify/validate a genome that is used for any
# comparison. The path to the source data, and a description of the genome are
# stored. We expect the path information to be live, so that a comparison may
# be run or re-run. The hash will be used to verify the contents of the file
# at the end of the path, when there is a run.
#
# Each pairwise comparison is stored (forward/reverse comparisons are stored
# separately, to allow method flexibility) in the comparisons table. The
# comparisons are tied directly to genomes, but only transitively to a
# particular run; this reduces redundancy, and allows pairwise comparison
# data to be used without recalculation, if the same input genome,
# path/hash, and analysis settings (fragsize for ANIb/ANIblastall and maxmatch
# for NUCmer),  are provided.
#
#    NOTE: Due to issues with Python/SQLite NULL queries,
#          fragsize default is zero
#
#
# The runs_genomes table provides a link so that each genome is associated with
# all runs in which it has participated, and each run can be associated with
# all the genomes that it has participated in.

# Create database tables
SQL_CREATEDB = """
   DROP TABLE IF EXISTS genomes;
   CREATE TABLE genomes (genome_id INTEGER PRIMARY KEY AUTOINCREMENT,
                         hash TEXT,
                         path TEXT,
                         length INTEGER,
                         description TEXT
                        );
   DROP TABLE IF EXISTS classes;
   CREATE TABLE classes (genome_id INTEGER NOT NULL,
                         run_id INTEGER NOT NULL,
                         class TEXT,
                         PRIMARY KEY (genome_id, run_id),
                         FOREIGN KEY(run_id) REFERENCES runs(run_id),
                         FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
                        );
   DROP TABLE IF EXISTS labels;
   CREATE TABLE labels (genome_id INTEGER NOT NULL,
                        run_id INTEGER NOT NULL,
                        label TEXT,
                        PRIMARY KEY (genome_id, run_id),
                        FOREIGN KEY(run_id) REFERENCES runs(run_id),
                        FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
                       );
   DROP TABLE IF EXISTS runs;
   CREATE TABLE runs (run_id INTEGER PRIMARY KEY AUTOINCREMENT,
                      method TEXT,
                      cmdline TEXT,
                      date TEXT,
                      status TEXT,
                      name TEXT
                     );
   DROP TABLE IF EXISTS runs_genomes;
   CREATE TABLE runs_genomes(run_id INTEGER NOT NULL,
                             genome_id INTEGER NOT NULL,
                             PRIMARY KEY (run_id, genome_id),
                             FOREIGN KEY(run_id) REFERENCES
                                                   runs(run_id),
                             FOREIGN KEY(genome_id) REFERENCES
                                                      genomes(genome_id)
                            );
   DROP TABLE IF EXISTS comparisons;
   CREATE TABLE comparisons (query_id INTEGER NOT NULL,
                             subject_id INTEGER NOT NULL,
                             aln_length INTEGER,
                             sim_errs INTEGER,
                             identity REAL,
                             cov_query REAL,
                             cov_subject REAL,
                             program TEXT,
                             version TEXT,
                             fragsize TEXT,
                             maxmatch TEXT,
                             PRIMARY KEY (query_id, subject_id, program,
                                          version, fragsize, maxmatch)
                            );
   DROP TABLE IF EXISTS runs_comparisons;
   CREATE TABLE runs_comparisons(run_id INTEGER NOT NULL,
                                 query_id INTEGER NOT NULL,
                                 subject_id INTEGER NOT NULL,
                                 program TEXT,
                                 version TEXT,
                                 fragsize TEXT,
                                 maxmatch TEXT,
                                 PRIMARY KEY (run_id, query_id, subject_id),
                                 FOREIGN KEY(run_id) REFERENCES
                                                       runs(run_id),
                                 FOREIGN KEY(query_id, subject_id, program,
                                             version, fragsize, maxmatch)
                                   REFERENCES
                                     comparisons(query_id, subject_id,
                                                 program, version, fragsize,
                                                 maxmatch)
                            );
   """

# Create indexes on hash in genome table
# The hash index is a standard single column index
# The hashpath index is a UNIQUE index to ensure that we don't duplicate the
# exact same file (though copies of the file in other places are allowed).
SQL_INDEXGENOMEHASH = """
   DROP INDEX IF EXISTS genomehash_index;
   CREATE INDEX genomehash_index ON genomes (hash);
   DROP INDEX IF EXISTS genomehashpath_index;
   CREATE UNIQUE INDEX genomehashpath_index ON genomes (hash, path);
"""

# Add an analysis run to the database
SQL_ADDRUN = """
   INSERT INTO runs (method, cmdline, date, status, name)
            VALUES (?, ?, ?, ?, ?);
"""

# Get information for a single analysis run
SQL_GETRUN = """
   SELECT * FROM runs WHERE run_id=?;
"""

# Get method for a single analysis run
SQL_GETMETHOD = """
   SELECT method FROM runs WHERE run_id=?;
"""

# Add a genome to the database
SQL_ADDGENOME = """
   INSERT INTO genomes (hash, path, length, description) VALUES (?, ?, ?, ?);
"""

# Associate a run with a genome
SQL_ADDRUNGENOME = """
   INSERT INTO runs_genomes (run_id, genome_id) VALUES (?, ?);
"""

# Get a specific genome hash
SQL_GETGENOMEHASH = """
   SELECT * FROM genomes WHERE hash=?;
"""

# Add a genome label to the database
SQL_ADDGENOMELABEL = """
   INSERT INTO labels (genome_id, run_id, label) VALUES (?, ?, ?);
"""

# Add a genome class to the database
SQL_ADDGENOMECLASS = """
   INSERT INTO classes (genome_id, run_id, class) VALUES (?, ?, ?);
"""

# Get a specific genome path
SQL_GETGENOMEPATH = """
   SELECT path FROM genomes WHERE genome_id=?;
"""

# Get the length of a genome
SQL_GETGENOMELENGTH = """
   SELECT length FROM genomes WHERE genome_id=?;
"""

# Get a specific genome hash/path combination
SQL_GETGENOMEHASHPATH = """
   SELECT * FROM genomes WHERE hash=? AND path=?;
"""

# Get all genome IDs associated with a specified run
SQL_GETGENOMESBYRUN = """
   SELECT genome_id FROM runs_genomes WHERE run_id=?;
"""

# Add a comparison to the database
SQL_ADDCOMPARISON = """
   INSERT INTO comparisons
     (query_id, subject_id, aln_length, sim_errs, identity, cov_query,
      cov_subject, program, version, fragsize, maxmatch)
     VALUES
     (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
"""

# Add a comparison to the database
SQL_ADDCOMPARISONLINK = """
   INSERT INTO runs_comparisons
     (run_id, query_id, subject_id, program, version, fragsize, maxmatch)
     VALUES (?, ?, ?, ?, ?, ?, ?);
"""

# Get a comparison (if it exists)
SQL_GETCOMPARISON = """
   SELECT * FROM comparisons WHERE query_id=? AND subject_id=? AND
                                   program=? AND version=? AND
                                   fragsize=? AND maxmatch=?;
"""

# Get all comparisons for a given run ID
SQL_GETRUNCOMPARISONS = """
   SELECT query.description, subject.description, comparisons.*
      FROM runs_comparisons
      JOIN comparisons USING
         (query_id, subject_id, program, version, fragsize, maxmatch)
      JOIN genomes query ON query_id=query.genome_id
      JOIN genomes subject ON subject_id=subject.genome_id
    WHERE runs_comparisons.run_id=?;
"""

# Get all analysis runs
SQL_GETALLRUNS = """
   SELECT run_id, name, method, date, cmdline FROM runs;
"""

# Get all genomes from the database
SQL_GETALLGENOMES = """
   SELECT genome_id, description, path, hash, length FROM genomes;
"""

# Get the JOIN of all runs to all genomes in the database
SQL_GETRUNSGENOMES = """
   SELECT runs.run_id, runs.name, runs.method, runs.date,
          genomes.genome_id, genomes.description, genomes.path,
          genomes.hash, genomes.length, classes.class, labels.label
     FROM runs
          JOIN runs_genomes ON runs.run_id=runs_genomes.run_id
          JOIN genomes ON runs_genomes.genome_id=genomes.genome_id
          JOIN classes ON runs.run_id=classes.run_id AND
                          genomes.genome_id=classes.genome_id
          JOIN labels ON runs.run_id=labels.run_id AND
                         genomes.genome_id=labels.genome_id
     ORDER BY runs.run_id, genomes.genome_id
"""

# Get the JOIN of all genomes to all runs in the database
SQL_GETGENOMESRUNS = """
   SELECT genomes.genome_id, genomes.description, genomes.path,
          genomes.hash, genomes.length,
          runs.run_id, runs.name, runs.method, runs.date
     FROM genomes, runs_genomes, runs
     WHERE genomes.genome_id=runs_genomes.genome_id AND
           runs_genomes.run_id=runs.run_id
     ORDER BY genomes.genome_id
"""

# Get the label for a genome in a given run
SQL_GETGENOMELABEL = """
   SELECT genomes.genome_id, genomes.hash, labels.label
     FROM genomes, labels
     WHERE genomes.genome_id=labels.genome_id AND
           labels.genome_id=? AND
           labels.run_id=?
     ORDER BY genomes.genome_id
"""

# Get the label for a genome in a given run
SQL_GETGENOMECLASS = """
   SELECT genomes.genome_id, genomes.hash, classes.class
     FROM genomes, classes
     WHERE genomes.genome_id=classes.genome_id AND
           classes.genome_id=? AND
           classes.run_id=?
     ORDER BY genomes.genome_id
"""

# Relabel a genome in a given run
SQL_UPDATEGENOMELABEL = """
   UPDATE labels
     SET label = ?
     WHERE genome_id = ? AND run_id = ?
"""

# Relabel a genome in a given run
SQL_UPDATEGENOMECLASS = """
   UPDATE classes
     SET class = ?
     WHERE genome_id = ? AND run_id = ?
"""


# DATABASE INTERACTIONS
# =====================


# Create an empty pyani SQLite3 database
def create_db(path):
    """Create an empty pyani SQLite3 database at the passed path."""
    conn = sqlite3.connect(path)
    with conn:
        cur = conn.cursor()
        cur.executescript(SQL_CREATEDB)
        cur.executescript(SQL_INDEXGENOMEHASH)


# Add a new run to the database
def add_run(dbpath, method, cmdline, date, status, name):
    """Add run information to the passed database, and return a run ID."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_ADDRUN, (method, cmdline, date, status, name))
    return cur.lastrowid


# Get a specified run row from the database
def get_run(dbpath, run_id):
    """Return the information associated with a specified run ID."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETRUN, (run_id,))
    return cur.fetchone()


# Get a specified run's method from the database
def get_method(dbpath, run_id):
    """Return the method used for a specified run."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETMETHOD, (run_id,))
    return cur.fetchone()[0]


# Add a new genome to the database
def add_genome(dbpath, genomehash, filepath, length, desc):
    """Add a genome to the passed SQLite3 database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        # The following line will fail if the genome is already in the
        # database, i.e. if the hash is not unique
        cur.execute(SQL_ADDGENOME, (genomehash, filepath, length, desc))
    return cur.lastrowid


# Associate a run ID with a genome ID
def add_genome_to_run(dbpath, run_id, genome_id):
    """Associate a run with a genome."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_ADDRUNGENOME, (run_id, genome_id))
    return cur.lastrowid


# Return the row corresponding to a single genome, defined by hash
def get_genome(dbpath, genomehash, path=None):
    """Return genome data if the passed hash is in the genomes table."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        if path is None:
            cur.execute(SQL_GETGENOMEHASH, (genomehash,))
        else:
            cur.execute(SQL_GETGENOMEHASHPATH, (genomehash, path))
        result = cur.fetchall()
    return result


# Return the filepath associated with a genome_id
def get_genome_path(dbpath, genome_id):
    """Return the file path associated with a genome_id."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMEPATH, (genome_id,))
        result = cur.fetchone()
    return result[0]


# Return the total length of a genome, identified by genome_id
def get_genome_length(dbpath, genome_id):
    """Return the genome length associated with a genome_id."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMELENGTH, (genome_id,))
        result = cur.fetchone()
    return result[0]


# Return genome IDs associated with a specific run
def get_genome_ids_by_run(dbpath, run_id):
    """Return list of genome IDs corresponding to the run with passed ID."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMESBYRUN, (run_id,))
        result = cur.fetchall()
    return [gid[0] for gid in result]


# Add a comparison to the database
def add_comparison(
    dbpath,
    qid,
    sid,
    aln_len,
    sim_errs,
    pid,
    qcov,
    scov,
    program,
    version,
    fragsize=0,
    maxmatch=None,
):
    """Add a single pairwise comparison to the database.

    NOTE: Due to issues with Python/SQLite NULL queries,
          fragsize default is zero
    """
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(
            SQL_ADDCOMPARISON,
            (
                qid,
                sid,
                aln_len,
                sim_errs,
                pid,
                qcov,
                scov,
                program,
                version,
                fragsize,
                maxmatch,
            ),
        )
    return cur.lastrowid


# Add a comparison/run link to the database
def add_comparison_link(
    dbpath, run_id, qid, sid, program, version, fragsize=0, maxmatch=None
):
    """Add a single pairwise comparison:run ID link to the database.

    NOTE: Due to issues with Python/SQLite NULL queries,
          fragsize default is zero
    """
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(
            SQL_ADDCOMPARISONLINK,
            (run_id, qid, sid, program, version, fragsize, maxmatch),
        )
    return cur.lastrowid


# Check if a comparison has been performed
def get_comparison(dbpath, qid, sid, program, version, fragsize=0, maxmatch=None):
    """Return the genome ID of a specified comparison.

    NOTE: Due to issues with Python/SQLite NULL queries,
          fragsize default is zero
    """
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETCOMPARISON, (qid, sid, program, version, fragsize, maxmatch))
        result = cur.fetchone()
    return result


# Get comparisons associated with a run
def get_comparisons_by_run(dbpath, run_id):
    """Return the comparisons associated with a specific run
    """
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETRUNCOMPARISONS, (run_id,))
        result = cur.fetchall()
    return result


# Return a table of all runs in the database
def get_all_runs(dbpath):
    """Return a table of all runs in the database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETALLRUNS)
        result = cur.fetchall()
    return result


# Return a table of all runs in the database
def get_all_genomes(dbpath):
    """Return a table of all genomes in the database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETALLGENOMES)
        result = cur.fetchall()
    return result


# Return a table of all runs in the database, joined with all genomes in the
# run
def get_genomes_by_runs(dbpath):
    """Return a table JOINing all runs to the genomes in the run."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETRUNSGENOMES)
        result = cur.fetchall()
    return result


# Return a table of all genomes in the database, joined with all the runs
# in which it participates
def get_runs_by_genomes(dbpath):
    """Return a table JOINing all genomes to the runs in the database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMESRUNS)
        result = cur.fetchall()
    return result


# Add a genome label to the database
def add_genome_label(dbpath, genome_id, run_id, label):
    """Add a single genome label to the database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_ADDGENOMELABEL, (int(genome_id), int(run_id), str(label)))
    return cur.lastrowid


# Relabel genome by genome ID and run ID
def update_genome_label(dbpath, genome_id, run_id, label):
    """Relabels a single genome in a given run."""

    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_UPDATEGENOMELABEL, (str(label), int(genome_id), int(run_id)))
    return cur.lastrowid


# Get a genome's label (for a specific run)
def get_genome_label(dbpath, genome_id, run_id):
    """Return the label associated with a genome for a given run."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMELABEL, (genome_id, run_id))
        result = cur.fetchone()
    if result is not None:
        return result[2]
    else:
        return "run {0}, genome {1}".format(run_id, genome_id)


# Add a genome class to the database
def add_genome_class(dbpath, genome_id, run_id, gclass):
    """Add a single genome class to the database."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_ADDGENOMECLASS, (int(genome_id), int(run_id), str(gclass)))
    return cur.lastrowid


# Reclass genome by genome ID and run ID
def update_genome_class(dbpath, genome_id, run_id, label):
    """Reclass a single genome in a given run."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_UPDATEGENOMECLASS, (str(label), int(genome_id), int(run_id)))
    return cur.lastrowid


# Get a genome's class (for a specific run)
def get_genome_class(dbpath, genome_id, run_id):
    """Return the class associated with a genome for a given run."""
    conn = sqlite3.connect(dbpath)
    with conn:
        cur = conn.cursor()
        cur.execute(SQL_GETGENOMECLASS, (genome_id, run_id))
        result = cur.fetchone()
    if result is not None:
        return result[2]
    else:
        return "run {0}, genome {1}".format(run_id, genome_id)


# RESULTS AS DATAFRAMES
# =====================

# Get comparisons for a single result, as a dataframe
def get_df_comparisons(dbpath, run_id):
    """Return a dataframe describing all comparisons for a single run."""
    data = get_comparisons_by_run(dbpath, run_id)
    headers = [
        "query",
        "subject",
        "query ID",
        "subject ID",
        "aligned length",
        "similarity errors",
        "percentage identity",
        "query coverage",
        "subject coverage",
        "program",
        "version",
        "fragsize",
        "maxmatch",
    ]
    df = pd.DataFrame(data)
    df.columns = headers
    return df


# Get all runs for which each genome is involved
def get_df_genome_runs(dbpath):
    """Return dataframe describing all runs for which a genome is involved."""
    data = get_runs_by_genomes(dbpath)
    headers = [
        "genome ID",
        "description",
        "path",
        "MD5 hash",
        "genome length",
        "run ID",
        "name",
        "method",
        "date run",
    ]
    df = pd.DataFrame(data)
    df.columns = headers
    return df


# Get all genomes used for each run
def get_df_run_genomes(dbpath):
    """Return dataframe describing all genomes involved with each run."""
    data = get_genomes_by_runs(dbpath)
    headers = [
        "run ID",
        "name",
        "method",
        "date run",
        "genome ID",
        "description",
        "path",
        "MD5 hash",
        "genome length",
        "class",
        "label",
    ]
    df = pd.DataFrame(data)
    df.columns = headers
    return df


# Get all genomes in the database
def get_df_genomes(dbpath):
    """Return dataframe describing all genomes in the database."""
    data = get_all_genomes(dbpath)
    headers = ["genome ID", "description", "path", "MD5 hash", "genome length"]
    df = pd.DataFrame(data)
    df.columns = headers
    return df


# Get all runs in the database
def get_df_runs(dbpath):
    """Return dataframe describing all runs in the database."""
    data = get_all_runs(dbpath)
    headers = ["run ID", "name", "method", "date run", "command-line"]
    df = pd.DataFrame(data)
    df.columns = headers
    return df


# Relabel genomes in the database
def relabel_genomes_from_file(dbpath, relabelfname, run_id, force=False):
    """Relabel genomes in the database using names in the passed file."""
    newlabels = parse_labelfile(relabelfname)
    genomes = get_df_genomes(dbpath)
    genomes.set_index("MD5 hash", inplace=True)
    for genomehash in newlabels:
        genome_id = genomes.loc[genomehash]["genome ID"]
        if force:
            add_genome_label(dbpath, genome_id, run_id, str(newlabels[genomehash]))
        else:
            update_genome_label(dbpath, genome_id, run_id, str(newlabels[genomehash]))


# Change classes of genomes in the database
def reclass_genomes_from_file(dbpath, reclassfname, run_id, force=False):
    """Reclass genomes in the database using names in the passed file."""
    newclasses = parse_labelfile(reclassfname)
    genomes = get_df_genomes(dbpath)
    genomes.set_index("MD5 hash", inplace=True)
    for genomehash in newclasses:
        genome_id = genomes.loc[genomehash]["genome ID"]
        if force:
            add_genome_class(dbpath, genome_id, run_id, str(newclasses[genomehash]))
        else:
            update_genome_class(dbpath, genome_id, run_id, str(newclasses[genomehash]))


# Parse labelfile for genomes in database
def parse_labelfile(fname):
    """Parse relabelling information for genomes in the database.

    The file should have tab-separated entries in the format
    <hash>\t<label>.

    Returns a dictionary keyed by genome hash, with value being the new label.
    """
    newlabels = {}
    with open(fname, "r") as infh:
        for line in infh:
            cleanline = line.strip().split("\t")
            newlabels[cleanline[0]] = cleanline[1]
    return newlabels


# RESULTS CLASS
# =============

# Class to produce/hold ANI results from a named run


class ANIResults(object):
    """Interfaces with the pyani database to extract and hold run output.

    Provides five dataframes, populated on instantiation

    - identity    %ID of aligned regions, symmetrical for ANIm,
                  not necessarily for other methods
    - coverage    %coverage of query sequence by aligned regions
    - aln_lengths total length of aligned region (symmetrical)
    - sim_errors  total number of alignment errors (mismatches, etc.
                  -symmetrical)
    - hadamard    dot product of coverage and identity matrices
    """

    def __init__(self, dbpath, run_id):
        """Initialise with empty dataframes and path to database."""
        self.dbpath = dbpath
        self.run_id = run_id
        self.method = get_method(self.dbpath, self.run_id)
        self.__get_labels()
        self.__initialise_dataframes()
        self.__get_data()

    def __initialise_dataframes(self):
        """Initialise empty dataframes for object.

        Creates five empty dataframes for ANI results.
        """
        self.identity = pd.DataFrame(
            index=self.genome_ids, columns=self.genome_ids, dtype=float
        )
        self.coverage = pd.DataFrame(
            index=self.genome_ids, columns=self.genome_ids, dtype=float
        )
        self.aln_lengths = pd.DataFrame(
            index=self.genome_ids, columns=self.genome_ids, dtype=float
        )
        self.sim_errors = pd.DataFrame(
            index=self.genome_ids, columns=self.genome_ids, dtype=float
        )
        self.hadamard = pd.DataFrame(
            index=self.genome_ids, columns=self.genome_ids, dtype=float
        )

    def __get_labels(self):
        """Retrieve genome IDs and labels for this run."""
        self.genome_ids = get_genome_ids_by_run(self.dbpath, self.run_id)
        self.labels = {
            genome_id: "_".join(
                [get_genome_label(self.dbpath, genome_id, self.run_id), str(genome_id)]
            )
            for genome_id in self.genome_ids
        }
        self.classes = {
            genome_id: "_".join(
                [get_genome_class(self.dbpath, genome_id, self.run_id), str(genome_id)]
            )
            for genome_id in self.genome_ids
        }
        self.lengths = {
            genome_id: get_genome_length(self.dbpath, genome_id)
            for genome_id in self.genome_ids
        }

    def __get_data(self):
        """Populate dataframes from database."""
        comparisons = get_df_comparisons(self.dbpath, self.run_id)
        for idx, row in comparisons.iterrows():
            qid, sid = row["query ID"], row["subject ID"]

            # Add percentage identity values; ANIm is symmetrical
            self.identity.loc[qid, sid] = row["percentage identity"]
            if self.method == "ANIm":
                self.identity.loc[sid, qid] = row["percentage identity"]

            # Add query coverage values; ANIm only has a single comparison,
            # so the lower triangle is subject coverage
            self.coverage.loc[qid, sid] = row["query coverage"]
            if self.method == "ANIm":
                self.coverage.loc[sid, qid] = row["subject coverage"]

            # Add similarity errors; ANIm is symmetrical
            self.sim_errors.loc[qid, sid] = row["similarity errors"]
            if self.method == "ANIm":
                self.sim_errors.loc[sid, qid] = row["similarity errors"]

            # Add alignment lengths; ANIm is symmetrical
            self.aln_lengths.loc[qid, sid] = row["aligned length"]
            if self.method == "ANIm":
                self.aln_lengths.loc[sid, qid] = row["aligned length"]

            # Calculate Hadamard matrix
            self.hadamard = self.identity * self.coverage

        # Populate remaining diagonals
        np.fill_diagonal(self.identity.values, 1.0)
        np.fill_diagonal(self.coverage.values, 1.0)
        np.fill_diagonal(self.sim_errors.values, 0.0)
        np.fill_diagonal(self.hadamard.values, 1.0)
        for idx, length in self.lengths.items():
            self.aln_lengths.loc[idx, idx] = length

        # Add labels and indices
        for matname in [
            "identity",
            "coverage",
            "aln_lengths",
            "sim_errors",
            "hadamard",
        ]:
            mat = getattr(self, matname)
            mat.index = [self.labels[val] for val in mat.index]
            mat.columns = [self.labels[val] for val in mat.columns]
