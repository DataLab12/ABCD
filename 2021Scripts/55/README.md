# status
Status implementation in C++
postprocessing implementation in python 

[statusMain.cpp](statusMainLinux.cpp) - main Data Lab edits for Linux here 
[statusMainWindows.cpp](statusMainWindows.cpp) - main Data Lab edits for Windows here

*Note on running any statusMain file*: Uncomment the one you wish to run (they are compiler dependent and comment out the others in your CMakeLists.txt (if using Cmake)

[data-test](https://git.txstate.edu/DataLab/data-test)  - test folder to test output 
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

* we will add more stats, one row per node
* consider adding ground truth column as first column
  * empty for nodes if no gt
  * original edge value for edge always 
  
[scripts] - py scripts for postprocessing output of statusMain.cpp 
  * Blane to add python files here 
  
[deprecated](deprecated) - old versions of the code

[status55.cpp](status55.cpp) - current ECL version 55  28 Aug 2020

