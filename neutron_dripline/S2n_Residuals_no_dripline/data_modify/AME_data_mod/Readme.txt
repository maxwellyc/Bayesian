To convert Atomic Mass Evaluation style outputs (.txt/.dat) into Excel sheets that ONLY contain data that's not from extrapolation, do the following:
THIS -- Python code AME_mod.py

1. Use Excel to open original outputs (.txt or .dat), choose the data file using ONLY Excel "Open" option, Text Import Wizard will ask you how do you want to read .dat/.txt files:
Step 1 - Choose Delimited
Step 2 - Select Tab AND Space
Click Finish

2. In Excel, copy all data (including label row) except the first column and paste into a new .dat/.txt file, because:
*** NOTE THAT THE 0 IN THE FIRST COLUMN AFTER MASS NUMBER 100 CAUSES SOME ISSUE BY MAKING THE FIRST COLUMN 0100, WHICH COMBINED WITH THE MASS NUMBER, YET TO BE FIXED. 04/10/18
***
   The first column of data is an annoying "0" or blank, it's purpose is to put a 0 in the first row of every new mass number A
   This makes it hard to read data properly using THIS to extract data properly since it uses "split", whenever there's a "0" instead of blank space, the array indices that correspond to data shifts back 1 position. 

3. Use the new .dat/.txt file as input for THIS, a file would generate that you can open in Excel neatly.
   
**Currently the indices "a*" (eg. as2n for 2-neutron separation) are set for S2n w/Error, S2p w/Error, to extract other data please change these indices.

Maxwell Y. Cao
2017.DEC.2  