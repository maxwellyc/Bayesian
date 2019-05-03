To create residue plots of 2-nucleon separation energies (and potentially other variables) theory vs experiment, do the following:

1. Construct input files for residue.py in ./Residue_plot/
   a. Modify AME data using AME_mod.py in ./AME_data_mod/ , follow instructions in 
      Readme.txt in that folder.
   b. Modify RMF data using RMF_S2_mod.py in ./RMF_data_mod/ , delete the header 
      paragraphs of the RMF files, leaving only the label row and the data, put
      the No-Header files in ./RMF_data_mod/RMFnoHeader/

2. Copy or move AME and RMF files into ./Residue_plot/data/, make sure you created corresponding arrays in residue.py when you add new masstables to plot 

3. Run residue.py