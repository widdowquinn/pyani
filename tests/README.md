# README.md - pyani/tests

This README describes the `pyani` testing suite

## Testing `pyani`

### Dependencies

The tests in this directory rely on the [`nose`](https://nose.readthedocs.org/en/latest/) package, which can be installed using

```
pip install nose
```

### Running tests

To run the tests in this directory with `nose` run the following command:

```
nosetests
```

This will run silently for quite a while (the comparisons are not quick), but should generate output that looks like this:

```

```


## The tests

### `test_concordance.py`

This tests the results of ANIm, ANIb, ANIblastall and TETRA analysis against the corresponding results from the [`JSpecies`](http://imedea.uib-csic.es/jspecies/) package.