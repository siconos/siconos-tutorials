This folder contains the code associated to the example studied in the paper: *Index-2 hybrid DAE: a case study with
well-posedness and numerical analysis*. The code simulating the example is separated into a C++ version and Python jupyter version.
Code executed and tested on Ubuntu 18.04.

## Requirements
 
- Both the python code and the C++ code require the installation of the Siconos toolbox: https://github.com/siconos/siconos

- The python code requires Jupyter notebook: https://jupyter.org/

- The Python code as been tested with a Python 3.6 core and require numpy, ipympl and matplotlib modules

- The C++ only require the same librairies as Siconos

## Python-Jupyter Code

This code as been used to generate the pictures in the papers. Please Refers to the comments in the notebook for correct usage.

## C++ code

This code behave in the same way as the Python code but is significantly faster.

- Use 

    ``siconos example_convergence.cpp`` 

to generate the results of the convergence tests (see usage in example_convergence.cpp comments).

- Use 

    ``siconos example_simulation.cpp``

to generate the results of a single simulation run.

