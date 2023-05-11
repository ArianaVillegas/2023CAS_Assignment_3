================================================================================
Complex Adaptive Systems Assignment 3
================================================================================

This repository holds the code for our (Ariana Villegas, Christopher Leap,
Emmanuel Ohiri) report for assignment 3.

Conda Environment
--------------------------------------------------------------------------------
Before running any of the code, create and activate our conda environment::

        $ conda env create -f=environment.yml
        $ conda activate antigenic


Mutation Accumulation
--------------------------------------------------------------------------------
To run the code for mutation accumulation, run::

        $ python -m mut_acc

This will perform a comparison between NextStrain database and the results 
of our method. All the figures generated in the report are stored in the 
figures folder inside mut_acc folder.
