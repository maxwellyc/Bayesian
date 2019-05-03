TO ADD NEW MASSTABLES:
1. Update masstable index list below AND in residue.py, currently:

   1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 
   7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,
   11-AME2003, 12-AME2016 


In residue.py:

2. Create file path f* at the beginning, eg.:

   f1 = open("data/SKMSall_nuclei.dat")

3. Create array l[*] to store f*.readline(), eg.:

   l[1] = f1.readlines()

4. Update index range i (theory), j (experiment) in for loops:
   a. Loops in section 2. Read in data from file
   b. Loops in section 3. Calculate residues
   c. Loops in section 4. Plotting:
      ResC loop, consists of S2n storage loop, S2p storage loop

5. Close file f* after section 3, eg.:
   
   f1.close()

6. Update masstable name file array mtN, currently:

   ["","SkM*","SkP","SLy4","SV-min","UNEDF0","UNEDF1","DD-ME2","DD-ME$\delta$",
     "DDPC1","NL3*","AME2003","AME2016"]
   
   Note that the first element is blank to match masstable index list

TO ADD NEW OBSERVABLES (XYZ) FOR PLOTTING:
1. Add additional elements or change elements in aME, aRMF, aAME that tells python which 
   column of the input masstable data file corresponds to XYZ
2. Create dictionary / arrays:
   
   Section 2: Add dictionary XYZ = {}, also in for loops here add XYZ[(N,Z,i)] = ...
   
   Section 3: a. Add dictionary ResXYZ = {},
              c. Create new LOOP BLOCK for XYZ residue calculation, follow example of S2n 
                 residue
              d. If you want to know difference in which nuclei is in masstable i but not 
                 in masstable j, follow example of newNuc, remeber to change k index of 
                 newNuc[(N,Z,k)], which is data type
   
   Section 4: Create new LOOP BLOCK for data type k of ResC[(Z,i,j,k,l)], follow example 
              of Data storage for S2n vs Neutron plot.

              k is data type, create copy of the first nested i,j,Z,N loop for new data 
              type;
              l index of ResC is x-axis data and y-axis data, later to be used by 
              numpy.plot, in S2n example, l=0 is Neutron number (has to be in increasing 
              order if you want to plot isotopic chain) to be plotted as x-position, l=1 
              is S2n residue to be plotted as y-position.
              If you want to change x,y axis data, modify:

              ResC[(Z,i,j,k,0)].extend([Y-AXIS VARIABLE])
              ResC[(Z,i,j,k,1)].extend([X-AXIS VARIABLE])

   Section 4.a.: Create ENTIRE new SUBSECTION like 4.a. S2N Plotting for plotting, mainly 
                 because you would need new axis labels, titles, fN, grids, ticks, texts,     
                 you DON'T necessary need: 
                 
                 if ( i==11 , j==12 ):
                 ... 
                 elif (i == 12, j == 11):
                 ...
                 
                 because these are for residue plots between the two experimental data 
                 sets, which needs different x/y axis range, y-labels and text location 
                 needs to change because x/y axis range changed

FOR BELOW OPTIONS, YOU NEED TO MANUALLY TUNE POSITION OF TICKS AND TEXT SUCH AS NUCLEI COUNT AND CHI-SQUARE

TO PLOT SELECT CHAINS OF ISOTOPE OR ISOTONE:
In section 4.0. Plotting options
Change values of zS, zE for S2n plots or nS, nE for S2p plots
You can, for example, plot Calcium chain of S2n residue data by specifying zS=20,zE=20, but the file name is still same as global plot, change it manually, i don't want to code for every potential option you'd want.

TO PLOT MASSTABLE 1 VS MASSTABLE 2, OR ADDITIONAL DATA TYPE k
or combination of such, change i,j,k in 4.0. Plotting options


MISCs:
1. fFormat: figure format type when saving
2. fN: file path and name when saving
3. pMod: change figure save directory between "test", "global","local", you are 
   responsible for creating actual directory since python won't do that for you.





