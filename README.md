# Enumeration of regular fractional factorial designs with four-level and two-level factors

This repository holds the code to reproduce the results from our paper "Enumeration of regular fractional factorial designs with four-level and two-level factors".

With this code you can run the three enumeration procedures mentioned in the paper:

- ST-NAUTY
- DOP-NAUTY
- MCS-Regular

To use the files in this project, first install the pipenv environement associated with it.
Assuming pipenv is installed on your machine (if follow [these instructions](https://pipenv-fork.readthedocs.io/en/latest/basics.html#)), run the following command:

```bash
pipenv install
pipenv shell
```

## Test results

Figure 1 of the paper display three graphs that illustrate the differences in running time between the three enumeration procedures.
To obtain the data used to generate the three graphs, run `test.sh`.
Then, to generate the three graphs, run the following command:

```bash
Rscript test_cases_individual_graph.R
```

The figures will be generated in the `figures/` folder.

## Usage example

The following example shows how to enumerate all 32-run designs of resolution III with one four-level factor and 6 two-level factors.

### Root design

You first need to initialize the *root* design, i.e. the design that contains the four-level factor and all the remaining basic two-level factors.
In this example, there are $log_2(32)-2*m=5-2=3$ remaining basic factors, so the root design is a $4^12^3$ design.
To generate the root design for a design with `32` runs, `1` four-level factor, with resolution `3` and using the search table (`ST`) method, run `init3.py` with the following options:

```bash
python init3.py -l 32 1 3 ST
```

### Extending the root design

To obtain a $4^12^6$ design, we still need to add three two-level factors to the root design.
To do this using the search table method, run the `extend3_st.py` script:

```bash
python extend3_st.py -l 32 1 4 3
python extend3_st.py -l 32 1 5 3
python extend3_st.py -l 32 1 6 3
```

or, for larger designs, use a `for` loop:

```bash
for i in {4..10}
do 
    python extend3_st.py -l 32 1 $i 3
done
```

The `-l` flag in the procedure means that logs are created, keeping track of the number of parents, candidates, and representatives.
