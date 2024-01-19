# ABCD: Algorithm For Balanced Component Discovery in Signed Networks

Authors: <em> Muhieddine Shebaro, Jelena Tešić </em>

![Highland Tribes Execution!](/images/animate.gif "Highland Tribes Clustering")

**graphC** (2024) is a scalable state-of-the-art hierarchical clustering algorithm for signed graphs capable of automatically detecting clusters without a predefined K hyperparameter (number of communities), no matrices, no decomposition spectral solvers, and no ground-truth labels. The algorithm is implemented in C++ and employs an efficient fundamental cycle basis discovery method to balance a connected component, performs Harary cuts, and selects the most optimal split for a connected component.

## Citation
Please cite the following publication: TBD

**BibTeX entry:**
```
TBD
```

## How to Run the Code 

* Simply download index.cpp and GraphBplus_Harary.cpp and compile the former source file like the following:

```
user:~$ g++ -fopenmp index.cpp
```
The signed graph must be in the following format to be compatible with graphC (src,dst,sign).
Preprocessing of the signed graph is embedded and neutral edges are treated as positive.

* To execute the compiled file, graphC utilizes 6 parameters in this order (iteration_count 𝛼 𝛽 ε time_limit γ):
```
user:~$ ./a.out input.txt 0.5 1 0.000001 -1 2
```
You can input -1 in time_limit to allow the algorithm to run the algorithm until it's finished. Minimum of γ is 2. Range of 𝛼 and 𝛽 is [0,1]. Range of ε is [0, infinity].
 
* graphC outputs 2 .txt files:

1. *_labels.txt: This contains the assigned clustering labels for each node (original node ID).
2. *_posin_negout.txt: This contains the history of change of the fraction of positive edges within communities and fraction of negative edges between communities including the overall improvement and unhappy ratio after each Harary split.

**Note:** Do <em> not</em> change the name of "GraphBplus_Harary.cpp" file. And if you run into errors related to the stack memory please run this command before executing the code:
```
user:~$ ulimit -s unlimited
```

[Data Lab @ TXST](DataLab12.github.io)


