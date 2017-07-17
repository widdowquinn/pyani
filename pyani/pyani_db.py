# -*- coding: utf-8 -*-
"""pyani_db.py

This module provides useful functions for creating and manipulating pyani's
SQLite3 databases

(c) The James Hutton Institute 2016-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2016-2017 The James Hutton Institute

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

import sqlite3

# SQL SCRIPTS
#==============
# The following is SQL for various database admin tasks, defined here to
# be out of the way when reading code.

# Create database tables
SQL_CREATEDB = """
   DROP TABLE IF EXISTS genomes;
   CREATE TABLE genomes (genome_id INTEGER PRIMARY KEY AUTOINCREMENT,
                         hash TEXT,
                         path TEXT,
                         name TEXT,
                         label TEXT,
                         class TEXT
                        );
   DROP TABLE IF EXISTS comparisons;
   CREATE TABLE comparisons (query_id INTEGER,
                             subject_id INTEGER,
                             date TEXT,
                             options TEXT,
                             description TEXT,
                             identity REAL,
                             coverage REAL,
                             mismatches REAL,
                             aligned_length REAL,
                             PRIMARY KEY (query_id, subject_id)
                            );
   """


# Create an empty pyani SQLite3 database
def create_db(path):
    """Create an empty pyani SQLite3 database at the passed path."""
    conn = sqlite3.connect(path)
    with conn:
        cur = conn.cursor()
        cur.executescript(SQL_CREATEDB)
