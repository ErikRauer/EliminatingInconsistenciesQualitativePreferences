# Supplementary Material for Eliminating Inconsistencies among CP-Theory Qualitative Preferences

This document contains a description of the various supplementary materials for the paper "Eliminating Inconsistencies among CP-Theory Qualitative Preferences".
These materials include the datasets, the experimental results, and the code for both generating the data sets and running the experiments.


## File Descriptions

- `minpreffbas_exact.py`: code used to implement and run the ILP based methods, including MFAS, MP, and MD problems
- `xmlfiletoprefgraphfileconverter.py`: code for reading and creating datasets, including reading xml files and converting to PFG usable by `minpreffbas_exact.py`, and converting CP-nets to CP-theories
- `CP-Nets`: folder containing dataset of CP-nets, sorted by number of variables, e.g. `n3d2g3` contains those with 3 variables. Contains xml files, PFG files, and output text.
- `CP-Theories`: folder containing dataset of CP-theories, sorted by number of variables, e.g. `n3d2g3` contains those with 3 variables. Contains xml files, PFG files, and output text.
- `CPNetExperimentalData.csv`: experimental results for the data set of CP-nets, see below for details
- `CPTheoryExperimentalData.csv`: experimental results for the data set of CP-theories, see below for details


## Running the Code

`minpreffbas_exact.py` requires the networkX and ortools libraries.
These can be installed by running
```
pip install networkx
```
and
```
pip install ortools
```


SL-IPGs are represented and read as text files where each line represents one edge.
The first two numbers are the tail and head vertices respectively and they're followed by a binary string representing which of the preference statements label the edge.

`minpreffbas_exact.py` iterates through files of the following form, running the MFAS, MD, MP using previous MD solution as a constraint, MP, and MD using previous MP solution as a constraint:
```
./CP-Nets/n{i}d2g3/n{i}d2g3-PFG{j}.txt
```
where 3 <= i <= 9 and 1<= j <= 3.
These are the text files corresponding with the CP-net data set and can be changed to run the CP-theory data set in a similar manner.
The corresponding output for each problem is printed to the files 
```
./CP-Nets/n{i}d2g3/n{i}d2g3-{j}-output.txt
```

This output first contains the data on the SL-IPG (number of preference statements, number of edges), then the solution data for the problems in the order described above.

The same is then done for the CP-theory data.


## Generating New Data

To generate new CP-nets, use https://cs.uky.edu/~goldsmit/papers/GeneratingCPnetCode.html which will generate them in a uniform manner.
The resulting preferences are in an XML format as described by https://www.ece.iastate.edu/~gsanthan/crisner.html
These can be converted to the format compatible with `minpreffbas_exact.py` via the functions 
```
xmlPrefFileToPGFconverter()
```
and
```
multiXmlPrefToPGFconverter()
```
in `xmlfiletoprefgraphfileconverter.py`.
Alternatively, one can first add relative importance variables to make them CP-theories using the function
```
CPNetToCPTheoryConverter()
```
also in `xmlfiletoprefgraphfileconverter.py`.


## Experimental Data Descriptions

The columns of `CPNetExperimentalData.csv` and `CPTheoryExperimentalData.csv` are described here.
The first group of rows each represent a CP-net or CP-theory depending on which file is being considered.
Unless otherwise noted, we will say CP-net below with the understanding that it is interchangeable with CP-theory when considering `CPTheoryExperimentalData.csv`.
The next set of rows are the average values grouped by number of preference variables, e.g. the average value for all CP-nets with 3 preference values.
The final row contains the total average across all CP-nets.
The columns are:

- \# Vars: number of variables in the CP-net
- \# Prefs: number of preferences in the CP-net
- \# Edges: number of edges in the CP-net
- MFAS/Min Prefs/Min Edges: indicater that the next set of columns are the results when running the ILP algorithms for the different problems. 
MFAS is the general objective of removing all cycles from the SL-IPG by minimizing the number of edges removed with no regard for preferences.
Min Prefs is the MP problem of removing all cycles from the SL-IPG and minimizing the number of preferences removed.
Min Edges is the MD problem of removing all cycles from the SL-IPG and minimizing the number of edges removed (while considering preferences in the removal).
The next 8-11 columns for all three problems are mostly the same, with the exception of 3 redundant or not applicable columns removed from MFAS:
    - Time (s): time that the run took in number of seconds
    - \# Prefs Rem: (Not in MFAS) number of preferences removed by the solution
    - \# Edges Rem: number of edges removed by the solution
    - \# Iters: number of iterations that the run took
    - Init \# Rep Cycs: initial number of representative cycles found for the first iteration
    - Final \# Rep Cycs: final number of representative cycles used in the final iteration when the solution was found
    - \# Diff: difference between Init \# Rep Cycs and Final \# Rep Cycs, how many representative cycles were added to the initial set of cycles
    - \% Diff: percentage that the set of representative cycles grew by, Final \# Rep Cycs divided by Init \# Rep Cycs
    - Cycs Per Iter: average number of cycles added each iteration, \# Diff divided by \# Iters
    - Init \# Constr: (Not in MFAS, redundant) initial number of constraints in the ILP solver used in the first iteration
    - Final \# Constr: (Not in MFAS, redundant) final number of constraints in the ILP solver for the final iteration when the solution was found
- The sequential objective problems are also presented in a similar manner as the normal problems above, except that the \# Edges or \# Preferences are split into two columns for the min and max of their ranges depending on the problem being considered
- Comparisons: Indicates that the remaining set of columns are for the comparison of different problem solutions
- MFAS w/: Comparisons between the number of edges removed by the MFAS solution and the other two problems:
    - \# Min Prefs: difference between the number of edges removed by the MP problem solution and the MFAS solution
    - \% Min Prefs: how many times more edges were removed by the MP problem solution than by the MFAS solution
    - \# Min Edges: difference between the number of edges removed by the MD problem solution and the MFAS solution
    - \% Min Edges: how many times more edges were removed by the MD problem solution than by the MFAS solution
- Min Pref w/ Min Edge: comparisons between the MP and MD solutions
    - \% Time (E/P): how many times longer the MD problem took than the MP problem
    - \# Prefs Rem (E-P): difference in number of preferences removed by the MD solution and the MP solution
    - \% Prefs Rem (E/P): how many times more preferences were removed by the MD solution than the MP solution
    - \# Edges Rem (P-E): (Not in CP-net data, since same as MFAS) difference in number of preferences removed by the MP solution and the MD solution
    - \% Edges Rem (E/P): (Not in CP-net data, since same as MFAS) how many times more edges were removed by the MP solution than the MD solution
    - \# Iterations (E-P): difference in number of iterations that the MD problem took compared to the MP problem
- A similar format is used for comparing Min Pref w/ MP then MD (denoted Both (P)), Min Edge w/ MP then MD, Min Pref w/ MD then MP (denoted Both (E)), and Min Edge w/ MD then MP. Note that all comparisons are made with respect to the minimum values of the ranges.
