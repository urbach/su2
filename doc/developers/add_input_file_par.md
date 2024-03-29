# Adding a parameter to the input file

In order to add a parameter that will be parsed from the input file, do the following steps.
For the details please look at the comments inside the source files. 

1. Decide to which ```struct``` it belongs to in ```su2/parameters.hh```. If it doesn't exist, create one for it.
2. Add a corresponding variable to the desired ```struct```.
3. If not already present, in ```su2/parse_input_file.hh``` add a declaration of the function parsing the input file for the desired program. The parsing functions have the same names but different namespaces according to the program they refer to.
4. In ```su2/parse_input_file.cc``` write or edit the definition of the corresponding parsing function.

The structure of the input file should be as hierachical as possible:
* resemblance with [tmLQCD](https://github.com/etmc/tmLQCD) 
* without global parameters dandling around
* without repetition of parameters (which may cause conflicts and undesired results)