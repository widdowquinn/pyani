.. _pyani-examples:

========
Examples
========

Using ``pyani index``
---------------------

``pyani index`` is used to created the required ``.md5`` files for each genome, as well as the ``classes.txt`` and ``labels.txt``, which will come in handy later on for labelling the datasets.  Say the directory ``./mygenomes/`` contains fasta-formated sequences of interest.  Running the following command will prepare the files for database construction

::
    pyani index ./mygenomes/ -v -l mygenomes_logfile.txt

Now, the directory will contain the original ``.fasta`` files as well as ``classes.txt`` and ``labels.txt``.  ``./mygenomes/` can now be used with the rest of the ``pyani`` workflow.



