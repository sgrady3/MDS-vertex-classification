Code to classify vertices in the context of the minimum dominating set (MDS) problem. The code is described in the paper "Domination-Based Classification Algorithms" that is pending publication. Three executables will be created, Classifier_A, Classifier_B, and Classiier_C. All three classifiers depend on an ILP solver. We utilized Gurobi (https://www.gurobi.com/) which is a world class solver and offers a free license for academic purposes. Classifier C also depends on a graph automorphism package. We utilized saucy (http://vlsicad.eecs.umich.edu/BK/SAUCY/) as it is highly tuned to sparse graphs. One can request the code for free for academic purposes. You will need to edit the Makefile to change include paths for your copy of Gurobi and Saucy. The methods are independent of ILP solver and graph automorphism finder, but the code will need to be edited to swap out these packages.
After exacutables have finished compiling, use via the command line. For example Classifier_A <graph-file.txt>. Output will be printed out to stdout. 

Description of algorithmic details can be found in the paper "Domination based classification algorithms for the controllability analysis of biological interaction networks" found here: https://www.nature.com/articles/s41598-022-15464-4
