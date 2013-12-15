cs294-davidson
==============

Jacobi-Davidson diagonalizer for full CI computation.

This should have the following parts:
1. Output file generated from Psi4, and a Python script that parses it into a 2-D (large) matrix
2. Some code to read the matrix in.
3. Some other code that diagonalizes the sucker, ala Jacobi-Davidson

About the python script
=======================
Run from terminal as

$> ./parse.py -i [outputfrompsi4].dat -b [#basis fcns]

Right now, the matrix is stored in memory, but I wasn't sure how to write it to a file.
