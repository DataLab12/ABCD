# status
Status implementation in C++ TO DO: 

* add more stats, one row per node
* consider adding ground truth column as first column
  * empty for nodes if no gt
  * original edge value for edge always 

[data-test](../data-test)  - test folder to test output 
  * Input: 
    * Raw_Data/*_edges.csv - test input files
  * Output: Output_Data/Features - created output files 
    * Output_Data/Features/*_nodeStat.csv 
'''
<tree no> NodeID, degree,status (normalized), authority (normalized),cumulative vacillation,
 E.g. 
1000 NodeID, degree, status,authority,vacillation,
0, 5. 0.65, 1.0, 2.63,
'''
	* Output_Data/Features/*_edgeStat.csv  
'''	
< tree no> From NodeID, To NodeID, span, vacillation (normalized),
E.g. 
< tree no> From NodeID, To NodeID, span, vacillation,
1000 0, 1, 0.35,0.58,
'''
  
[scripts] - py scripts for postprocessing output of statusMain.cpp 
	
##[55](55/) - last version we worked on
* [statusMain.cpp](55/statusMainLinux.cpp) - main Data Lab edits for Linux here 
* [statusMainWindows.cpp](55/statusMainWindows.cpp) - main Data Lab edits for Windows here
*  Note on running any statusMain file: Uncomment the one you wish to run (they are compiler dependent and comment out the others in your CMakeLists.txt (if using Cmake)


#graphB+ v0.2

* graphB+ is a balancing algorithm for signed social network graphs. It operates on graphs stored in CSV format. A sample input graph can be found [here](graph.csv).

* The inputs are simple text files consisting of a header line followed by an arbitrary number of lines, each specifying one edge of the graph, where the edge weight is the sign (either 1 or -1):

```From Node ID, To Node ID, Edge Weight   10,33,-1   8,5,1'''

* [graphBplus\_02.cpp](graphBplus_02.cpp) is the OpenMP C++ code (with g++ intrinsics). Note that graphB+ is protected by the 3-clause BSD license included in the beginning of the code.

* The OpenMP C++ code can be compiled as follows:
```g++ -O3 -march=native -fopenmp graphBplus_02.cpp -o graphBplus'''
* To run the code on the input file `graph.csv` with 100 samples and save the results of the balanced solutions in out.csv, enter:
```./graphBplus graph.csv 100 out.csv'''

* To obtain the inputs used in the paper listed below and convert them to our format, download and run the file [input\_script.sh](input_script.sh). Note that this script takes about an hour to execute and requires a large amount of disk space.

### Publication

G. Alabandi, J. Tesic, L. Rusnak, and M. Burtscher. "Discovering and Balancing Fundamental Cycles in Large Signed Graphs." Proceedings of the 2021 ACM/IEEE International Conference for High-Performance Computing, Networking, Storage and Analysis. November 2021. This work has been supported in part by the hardware donation from NVIDIA Corporation.



  
