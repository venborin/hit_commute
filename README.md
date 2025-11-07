# hit_commute
A small program to calculate the hitting and commute times as described in the paper of Chakra Chennubhotla and Ivet Bahar:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030172

Installation:
g++ -O3 -std=c++17 hit_commute.cpp -o hit_commute -I /path/to/eigen3/include

Usage:
./hit_commute  --pdb FILE.pdb                               # to get all pairwise hitting and commute time
./hit_commute  --pdb FILE.pdb --cutoff [FLOAT]              # same within a cutoff of [float]
./hit_commute  --pdb FILE.pdb --res CHAIN_ID:RESIDUE_NUMBER  # to get the HCs from residue CHAIN_ID:RESIDUE_NUMBER to all others
