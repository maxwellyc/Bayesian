###################################################################################################
### This Code is written to extract and plot residues of 2 nucleon separation energy from       ###
### theoretical masstable calculations with DFT vs experimental data                            ###
### EDFs includes:                                                                              ###
### Skyrme: SkM*, SkP, SLy4, SV-min, UNEDF0, UNEDF1                                             ###
### Skyrme data can be found on http://massexplorer.frib.msu.edu/content/masstables             ###
### Relativistic Mean-Field (RMF): DD-ME2, DD-MED, DDPC1, NL3*                                  ###
### RMF data refer to studies in: Physical Review C 89, 054320 (2014)                           ###
### Atomic Mass Evaluation (AME): AME2003, AME2016                                              ###
### AME data from AMDC: https://www-nds.iaea.org/amdc/                                          ###
### Maxwell Y. Cao 2017.Dec.2                                                                   ###
###################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
import scipy.interpolate
import scipy.ndimage as ndimage
import itertools as it
import math
import matplotlib.gridspec as gridspec

# For Adobe illustrator text
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

#For calculation of odd values please refer to: http://massexplorer.frib.msu.edu/content/masstables/Odd_Values_from_Even_Data.pdf
def residues():
 # function that detects if arg. is a number
 def isNum(s):
  try:
    float(s)
    return True
  except ValueError:
    return False
 #Atomic mass unit from Wikipedia
 # Note that FRDM2012 uses 1u = 931.5014MeV
 uAtomic = 931.49        # MeV rounded to 2 decimals
 #uAtomic = 931.4940954  #MeV
 #Proton and neutron masses From AME 2012, which HFB-24 was fitted to
 #M_P = 1.00782503223 * uAtomic  #938.27
 #M_N = 1.00866491585 * uAtomic  #939.57
 M_P = 938.27  # MeV rounded to 2 decimals
 M_N = 939.57  # MeV rounded to 2 decimals
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 17-TRIUMF/RIKEN                                                           !!#
 #!! 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2                                           !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #mtN stands for masstableName
 mtN = ["","SkM*","SkP","SLy4","SV-min","UNEDF0","UNEDF1","DD-ME2","DD-ME$\delta$","DDPC1","NL3*",
        "AME2003","AME2016","JYFLTRAP2017","FRDM2012","HFB24","UNEDF2","TRIUMF"]
 ###################################################################################################
 ###                                      CAUTION! WARNING!                                      ###
 ###                                      CAUTION! WARNING!                                      ###
 ###                                      CAUTION! WARNING!                                      ###
 ### ALL MASSTABLES CONTAINS ODD DATA, ODD DATA IN Relativistic Mean-Field (RMF) MASSTABLES      ###
 ### ARE EXTRAPOLATED WITH INCORRECT PAIRING GAPS, DO NOT USE ODD DATA FROM RMF MASSTABLES!!!    ###
 ### 2017.DEC.2                                                                                  ###
 ###################################################################################################
 
 ###################################################################################################
 ###                                    1. DEFINE FILE PATHS                                     ###
 ###################################################################################################
 #
 #mass explorer masstables
 dir = "data/max_raw_no_dripline/"
 f1 = open(dir+"SKMS_all_nuclei_max.dat"); f2 = open(dir+"SKP_all_nuclei_max.dat")
 f3 = open(dir+"SLY4_all_nuclei_max.dat"); f4 = open(dir+"SV-MIN_all_nuclei_max.dat")
 f5 = open(dir+"UNEDF0_all_nuclei_max.dat"); f6 = open(dir+"UNEDF1_all_nuclei_max.dat")
 f16 = open(dir+"UNEDF2_all_nuclei_max.dat")
 #column indices of Z, N, S2n, S2p, S1n, S1p. ME stands for MassExplorer
 aME = [1,2,8,6,7,5]
 #
 #RMF masstables, !!!DO NOT USE ODD DATA!!! CONTAINED IN RMF masstables
 #THESE ODD DATA ARE CALCULATED WITH WRONG PAIRING GAPS, ONLY USE EVEN S2N, S2P DATA
 f7 = open("data/ddme2-sep.dat"); f8 = open("data/ddmed-sep.dat")
 f9 = open("data/ddpc1-sep.dat"); f10 = open("data/nl3s-sep.dat")
 #column indices of Z, N, S2n, S2p. RMF stands for Relativistic Mean Field
 aRMF = [0,1,6,4]
 #additional data, as requested by collaborators:
 f14 = open("data/FRDM2012.dat"); f15 = open("data/hfb24.dat")
 #colum indices of Z, N, Binding energy for FRDM-2012
 aFRDM = [0,1,13]
 #colum indices of Z, A, Mass excess, S1n, S1p for HFB-24
 aHFB24 = [0,1,9,6,7]
 
 #2003 & 2016 AME masstables, S2n, S2p
 f11 = open("data/2003AME_S2.dat"); f12 = open("data/2016AME_S2.dat")
 #2003 & 2016 AME masstables, S1n, S1p
 f17= open("data/2003AME_S1.dat"); f18 = open("data/2016AME_S1.dat")
 #AME stands for Atomic Mass Evaluation, this aAME is good for both *_S2 file and *_S1 file
 #column indices of Z, N, S2n/S1n, S2n/S1n errors, S2p/S1p, S2p/S1p errors (error data only exists for AME)
 aAME = [0,1,2,3,4,5]
 #JYFLTRAP 2017 data
 f13 = open("data/rare_earthS2n.dat")
 #column indices of Z, N, mass excess, mass excess error/,S2n, S2n errors
 aJYFL = [0,1,2,3]#,4,5]
 
 l = {}
 # Skyrme masstables
 l[1] = f1.readlines(); l[2] = f2.readlines(); l[3] = f3.readlines(); l[4] = f4.readlines()
 l[5] = f5.readlines(); l[6] = f6.readlines();  l[16] = f16.readlines()
 # RMF masstables
 l[7] = f7.readlines(); l[8] = f8.readlines(); l[9] = f9.readlines(); l[10] = f10.readlines();
 # S2n,S2p of AME2003, AME2016; JYFLTRAP
 l[11] = f11.readlines(); l[12] = f12.readlines(); l[13] = f13.readlines();
 # FRDM-2012 & HFB-24
 l[14] = f14.readlines(); l[15] = f15.readlines();
 # S1n,S1p from AME2003, AME2016
 l[17] = f17.readlines(); l[18] = f18.readlines()
 
 
 
 ###################################################################################################
 ### S2* dictionaries contain all the separation energies read from the 12 masstables above      ###
 ### The key is structured as: (N,Z,masstable index),                                            ###
 ### Error dictionaries only contain errors from AME2003 and AME2016                             ###
 ### To aviod confusion, the masstable index of error dictionaries for AME2003/AME2016           ###
 ### is still 11/12, in the future it's possible to include errors for theories                  ###
 ### (Neutron #, Proton #, theory masstable index(1~10), exp masstable index (11~12)  )          ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN  !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 
 ###################################################################################################
 ###                                  2. READ IN DATA FROM FILE                                  ###
 ###################################################################################################
 S2n = {}; S2nErr = {}; S2p = {}; S2pErr = {}; BE={}; nMax = 2; zMax = 2 ; nMax1 = 2; zMax1 = 2
 S1n = {}; S1nErr = {}; S1p = {}; S1pErr = {}
 #low_lim is a filter for separation value, eg. when low_lim is 0 only non-negative ( >= 0) values will be recorded
 low_lim = -10000.0
 #Mass Explorer data read and store, use aME:
 for i in it.chain(range(1,7), range(16,17)):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aME[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aME[0])])+0.0001)      #Proton Number
     if (nMax1<N):
      nMax1 = N
     if (zMax1<Z):
      zMax1 = Z
     #S2n
     data = str(ss[int(aME[2])])
     if (isNum(data) and float(data) >= low_lim ):
      S2n[(N,Z,i)] = float(data)
     #S2p
     data = str(ss[int(aME[3])])
     if (isNum(data) and float(data) >= low_lim ):
      S2p[(N,Z,i)] = float(data)
     #S1n
     data = str(ss[int(aME[4])])
     if (isNum(data) and float(data) >= low_lim ):
      S1n[(N,Z,i)] = float(data)
     #S1p
     data = str(ss[int(aME[5])])
     if (isNum(data) and float(data) >= low_lim ):
      S1p[(N,Z,i)] = float(data)
    except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue
 #RMF data read and store, use aRMF:
 #There is no S1n S1p for RMF data, pairing gap not available.
 for i in range(7,11):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aRMF[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aRMF[0])])+0.0001)      #Proton Number
     if (nMax1<N):
      nMax1 = N
     if (zMax1<Z):
      zMax1 = Z
     #S2n
     data = str(ss[int(aRMF[2])])
     if (isNum(data)):
      S2n[(N,Z,i)] = float(data)
     #S2p
     data = str(ss[int(aRMF[3])])
     if (isNum(data)):
      S2p[(N,Z,i)] = float(data)
    except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue
 #Experimental AME data read and store, use aAME:
 #S2n & S2p:
 for i in range(11,13):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aAME[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aAME[0])])+0.0001)      #Proton Number
     #S2n
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     data = str(ss[int(aAME[2])])
     if (isNum(data) and float(data) >= 0.0 ):
      S2n[(N,Z,i)] = float(data)#/1000.0             #If experimental data is in KeV, convert to MeV
      #S2n Error, there's an indent here, only store error if there's S2n data
      data = str(ss[int(aAME[3])])
      if (isNum(data)):
       S2nErr[(N,Z,i)] = float(data)#/1000.0         #If experimental data is in KeV, convert to MeV
     #S2p
     data = str(ss[int(aAME[4])])
     if (isNum(data) and float(data) >= 0.0 ):
      S2p[(N,Z,i)] = float(data)#/1000.0             #If experimental data is in KeV, convert to MeV
      #S2p Error, there's an indent here, only store error if there's S2p data
      data = str(ss[int(aAME[5])])
      if (isNum(data)):
       S2pErr[(N,Z,i)] = float(data)#/1000.0         #If experimental data is in KeV, convert to MeV
    except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue
 #S1n & S1p:
 for i in range(17,19):
  for line in l[i]:
    # This i to ii ensure that S1n[(N,Z,11)] is from AME2003, S1n[(N,Z,12)] from AME2016, to avoid confusion
    # The difference in index for experimental S1n is contained in the data reading stage only
    ii = 0
    if (i == 17) : ii = 11
    if (i == 18) : ii = 12
    ss = line.split()
    try:
     N = int(float(ss[int(aAME[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aAME[0])])+0.0001)      #Proton Number
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     data = str(ss[int(aAME[2])])
     if (isNum(data) and float(data) >= 0.0 ):
      S1n[(N,Z,ii)] = float(data)#/1000.0             #If experimental data is in KeV, convert to MeV
      #S1n Error, there's an indent here, only store error if there's S1n data
      data = str(ss[int(aAME[3])])
      if (isNum(data)):
       S1nErr[(N,Z,ii)] = float(data)#/1000.0         #If experimental data is in KeV, convert to MeV
     data = str(ss[int(aAME[4])])
     if (isNum(data) and float(data) >= 0.0 ):
      S1p[(N,Z,ii)] = float(data)#/1000.0             #If experimental data is in KeV, convert to MeV
      #S1n Error, there's an indent here, only store error if there's S1n data
      data = str(ss[int(aAME[5])])
      if (isNum(data)):
       S1pErr[(N,Z,ii)] = float(data)#/1000.0         #If experimental data is in KeV, convert to MeV
    except (ValueError, IndexError):        #N,Z, or, S1n are not numbers
       continue
 #JYFLTRAP 2017
 for i in range(13,14):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aJYFL[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aJYFL[0])])+0.0001)      #Proton Number
     #S2n
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     #S2n and error
     data = str(ss[int(aJYFL[2])])
     if (isNum(data)):
      S2n[(N,Z,i)] = float(data)/1000.0             #Experimental data is in KeV, converting to MeV
      data = str(ss[int(aJYFL[3])])
      if (isNum(data)):
       S2nErr[(N,Z,i)] = float(data)/1000.0         #Experimental data is in KeV, converting to MeV
     #S1n and error
     data = str(ss[int(aJYFL[4])])
     if (isNum(data)):
      S1n[(N,Z,i)] = float(data)/1000.0             #Experimental data is in KeV, converting to MeV
      data = str(ss[int(aJYFL[5])])
      if (isNum(data)):
       S1nErr[(N,Z,i)] = float(data)/1000.0         #Experimental data is in KeV, converting to MeV
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 #FRDM 2012, USE aFRDM
 for i in range(14,15):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aFRDM[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aFRDM[0])])+0.0001)      #Proton Number
     #S2n
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     data = str(ss[int(aFRDM[2])])
     if (isNum(data)):
      BE[(N,Z,i)] = -1.0 * float(data)              #Binding energy
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 #HFB-24, USE aHFB24
 for i in range(15,16):
  for line in l[i]:
    ss = line.split()
    try:
     Z = int(float(ss[int(aHFB24[0])]))      #Neutron Number
     A = int(float(ss[int(aHFB24[1])]))      #Mass Number
     N = A - Z
     #S2n
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     data = str(ss[int(aHFB24[2])])                 #For HFB-24, input is mass excess
     if (isNum(data)):
      #Calculate Binding Energy from mass excess
      BE[(N,Z,i)] = round(-1.0 * ( Z*M_P + N*M_N - ( float(data) + A*uAtomic ) ), 2)
     data = str(ss[int(aHFB24[3])])
     # In HFB-24 masstable, some values are marked 999.99, I believe these are invalid
     if (isNum(data)) and data != "999.99":
      if (float(data) >= low_lim ):
       S1n[(N,Z,i)] = round(float(data),2)
     data = str(ss[int(aHFB24[4])])
     # In HFB-24 masstable, some values are marked 999.99, I believe these are invalid
     if (isNum(data)) and data != "999.99":
      if (float(data) >= low_lim ):
       S1p[(N,Z,i)] = round(float(data),2)
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 # Manually adding TRIUMF 54,55,56 titanium data points:
 # https:// journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.062503
 S2n[(32,22,17)] = 12.40794 ; S2nErr[(32,22,17)] = 0.01628
 S1n[(33,22,17)] = 4.159318 ; S1nErr[(33,22,17)] = 0.033121
 S2n[(34,22,17)] = 9.7187   ; S2nErr[(34,22,17)] = 0.1223
 S1p[(,,17)] = ; S1pErr[(,,17)] =
 S2p[(,,17)] = ; S2pErr[(,,17)] =
 # Manually adding 55,56,57 calcium masses calculated S1/2n data points, RIKEN
 #https:// journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022506
 S1n[(35,20,17)] = 1.560733 ; S1nErr[(35,20,17)] = 0.167171
 S2n[(36,20,17)] = 4.492051 ; S2nErr[(36,20,17)] = 0.254649
 S1n[(37,20,17)] = 1.931318 ; S1nErr[(37,20,17)] = 1.021078
 S1p[(,,17)] = ; S1pErr[(,,17)] =
 S2p[(,,17)] = ; S2pErr[(,,17)] =
 #
 ###################################################################################################
 ###          2.b. CALCULATE SEPARATION ENERGY, FOR INPUT DATA ONLY CONTAINING MASSES            ###
 ###################################################################################################
 # HFB-24 and FDRM2012 only contain mass excess / binding energy data, calculating separation energies:
 # FRDM-2012 & HFB-24
 d_cache = 0.0
 for i in range(14,16):
  for Z in range(2,zMax+1):
   for N in range(2,nMax+1):
    if ( ( (N,Z,i) in BE) and ( (N-2,Z,i) in BE) ):
     d_cache =  - BE[(N,Z,i)] + BE[(N-2,Z,i)]
     if d_cache >= low_lim : S2n[(N,Z,i)] = d_cache
    if ( ( (N,Z,i) in BE) and ( (N,Z-2,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-2,i)]
     if d_cache >= low_lim : S2p[(N,Z,i)] = d_cache
 # Only FRDM-2012 needs to calculate S1n/S1p from Binding
 for i in range(14,15):
  for Z in range(2,zMax+1):
   for N in range(2,nMax+1):
    if ( ( (N,Z,i) in BE) and ( (N-1,Z,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N-1,Z,i)]
     if  d_cache >= low_lim : S1n[(N,Z,i)] = d_cache
    if ( ( (N,Z,i) in BE) and ( (N,Z-1,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-1,i)]
     if  d_cache >= low_lim : S1p[(N,Z,i)] = d_cache

 ###################################################################################################
 ###                    3. CALCULATE RESIDUES, RESULTS STORED IN ResS2n, ResS2p                  ###
 ###################################################################################################

 ###################################################################################################
 ### ResS2n, ResS2p are Dictionaries for S2n, S2p residues between theory and exp.               ###
 ### The key is structrued as (N,Z,i,j):                                                         ###
 ### (Neutron #, Proton #, theory masstable index(1~10), exp masstable index (11~12)  )          ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 ResS2n = {}; ResS2p = {}; ResS1n = {}
 rms_s2n = {}; rms_s1n = {}
 zs = 2; ns = 2
 res1_p = {}; res2_p = {}
 del_switch = 0
 # delete from del_list1 is what we should've done in paper 1's RMS
 del_list1 = {(98,82,1),(104,84,2),(104,84,5),(8,2,11),(16,8,11),(32,20,11),(64,38,11),(84,50,11),(22,12,11)}

 # delete from del_list1 to produce paper 1's RMS
 del_list2 = {(98,82,1),(104,84,2),(104,84,5),(8,2,11),(16,8,11),(32,20,11),(22,12,11)}

 if del_switch == 1:
  for d in del_list1:
   if d in S2n:
    del S2n[d]
#  S2n[(98 ,60,12)] = S2n[(98 ,60,13)]
#  S2n[(100,62,12)] = S2n[(100,62,13)]
#  S2n[(100,64,12)] = S2n[(100,64,13)]
#  S2n[(102,64,12)] = S2n[(102,64,13)]
 elif del_switch == 2:
  for d in del_list2:
   if d in S2n:
    del S2n[d]

 
 # model 1~10, 14,15,16; exp 11~13
 for i in it.chain(range(1,7), range(14,17)):
  res2_p[mtN[i]] = 0
  for j in it.chain(range(11,14),[17]):
    for Z in range(2,zMax+1):
     for N in range(2,nMax+1):
      #S2n residue calculation
      if ( (N,Z,i) in S2n and (N,Z,j) in S2n ):
       ResS2n[(N,Z,i,j)] = S2n[(N,Z,j)] - S2n[(N,Z,i)]
       if j == 11 and (not N%2) and (not Z%2): res2_p[mtN[i]] = res2_p[mtN[i]] + 1
       #S2p residue calculation
      if ( (N,Z,i) in S2p and (N,Z,j) in S2p ):
       ResS2p[(N,Z,i,j)] = S2p[(N,Z,j)] - S2p[(N,Z,i)]

 # model 1~6, 14,15,16; exp 11~13
 for i in it.chain(range(1,7), range(14,17)):
#  res1_p[mtN[i]] = 0
  for j in it.chain(range(11,14),[17]):
    for Z in range(2,zMax+1):
     for N in range(2,nMax+1):
      #S1n residue calculation
      if ( (N,Z,i) in S1n and (N,Z,j) in S1n ):
       ResS1n[(N,Z,i,j)] = S1n[(N,Z,j)] - S1n[(N,Z,i)]

 # Find rms of S2n points in AME2016 but not in AME2003, starting from zs, ns
 
 #Only Skyrme and FRDM-2012, HFB-24 has S1n values.
 for i in it.chain(range(1,7), range(14,17)):
  print (mtN[i])
  rms_s2n[(i,12)] = 0
  rms_s1n[(i,12)] = 0
  count_rms = 0
  count_rms1= 0
  #S2n
  for Z in range(zs,zMax+1,2):
   for N in range(ns,nMax+1,2):
    if (N,Z,i,12) in  ResS2n and (N,Z,i,11) not in ResS2n :
     rms_s2n[(i,12)] = rms_s2n[(i,12)] + ResS2n[(N,Z,i,12)]**2
     count_rms = count_rms + 1
     #print(Z,N,round(ResS2n[(N,Z,i,12)],6))
  rms_s2n[(i,12)] = math.sqrt(rms_s2n[(i,12)] / float(count_rms))

  print (count_rms)
  print (rms_s2n[(i,12)])
  print ("**************")

  #S1n
  for Z in range(zs,zMax+1):
   for N in range(ns,nMax+1):
    if (N,Z,i,12) in  ResS1n and (N,Z,i,11) not in ResS1n :
     rms_s1n[(i,12)] = rms_s1n[(i,12)] + ResS1n[(N,Z,i,12)]**2
     count_rms1 = count_rms1 + 1
  rms_s1n[(i,12)] = math.sqrt(rms_s1n[(i,12)] / float(count_rms1))
 
 # rms of AME2003 with Z>=zs, N>=zs
 #S2n
 for i in it.chain(range(1,7), range(14,17)):
  rms_s2n[(i,11)] = 0
  count_rms = 0
  for Z in range(zs,zMax+1,2):
   for N in range(ns,nMax+1,2):
    if (N,Z,i,11) in  ResS2n:
     rms_s2n[(i,11)] = rms_s2n[(i,11)] + ResS2n[(N,Z,i,11)]**2
     count_rms = count_rms + 1
  rms_s2n[(i,11)] = math.sqrt(rms_s2n[(i,11)] / float(count_rms))
#  print (mtN[i], '/ 2003')
#  print (rms_s2n[(i,11)])


 #S1n
 for i in it.chain(range(1,7), range(14,17)):
  rms_s1n[(i,11)] = 0
  count_rms1= 0
  for Z in range(zs,zMax+1):
   for N in range(ns,nMax+1):
    if (N,Z,i,11) in  ResS1n:
     rms_s1n[(i,11)] = rms_s1n[(i,11)] + ResS1n[(N,Z,i,11)]**2
     count_rms1 = count_rms1 + 1

 # newNucS2n records (N,Z) that has S2n or S2p values in AME2016 but not in AME2003, key is (N,Z,k)
 # k is data type, k=0: S2n, k=1: S2p
 newNucS2n = {}; newNucS1n = {}
 for i in range(12,14):
  for j in range(11,13):
    #S2n new in 2016 but not in 2003
    for Z in range(2,zMax+1,2):
     for N in range(2,nMax+1,2):
      #Records new nuclei with S2n data in 2016 compared to 2013
      if (i == 12 and j == 11):
       if ( (N,Z,i) in S2n and (N,Z,j) not in S2n ):
        newNucS2n[(N,Z,0)] = 1
       elif ( (N,Z,i) in S2n ):
        newNucS2n[(N,Z,0)] = -1
      #Records new nuclei with S2p data in 2016 compared to 2013
      if (i == 12 and j == 11):
       if ((N,Z,i) in S2p and (N,Z,j) not in S2p):
        newNucS2n[(N,Z,1)] = 1
       elif ( (N,Z,i) in S2p ):
        newNucS2n[(N,Z,1)] = -1
      #Records new nuclei with S2n data in 2017 compared to 2016:
      if (i == 13 and j == 12 ):
       if ( (N,Z,i) in S2n and (N,Z,j) not in S2n ):
        newNucS2n[(N,Z,0)] = 2
      if (i == 13 and j == 12 ):
       if ( (N,Z,i) in S2p and (N,Z,j) not in S2p ):
        newNucS2n[(N,Z,1)] = 2
    #S1n new in 2016 but not in 2003
    for Z in range(2,zMax+1):
     for N in range(2,nMax+1):
      #Records new nuclei with S1n data in 2016 compared to 2013
      if (i == 12 and j == 11):
       if ( (N,Z,i) in S1n and (N,Z,j) not in S1n ):
        newNucS1n[(N,Z,0)] = 1
       elif ( (N,Z,i) in S1n ):
        newNucS1n[(N,Z,0)] = -1
      #Records new nuclei with S1n data in 2017 compared to 2016:
      if (i == 13 and j == 12 ):
       if ( (N,Z,i) in S1n and (N,Z,j) not in S1n ):
        newNucS1n[(N,Z,0)] = 2


 f1.close(); f2.close(); f3.close(); f4.close(); f5.close(); f6.close(); f7.close(); f8.close()
 f9.close(); f10.close(); f11.close(); f12.close(); f13.close(); f14.close(); f15.close()
 f16.close(); f17.close(); f18.close()

 ###################################################################################################
 ###                                         4. PLOTTING                                         ###
 ###################################################################################################

 ###################################################################################################
 ### ResC: Dictionary of arrays for S2n plotting, the key is structured as: (Z,i,j,k,l):         ###
 ### ResC: Dictionary of arrays for S2n plotting, the key is structured as: (Z,i,j,k,l):         ###
 ### Z: proton #;                                                                                ###
 ### i: theory masstable index 1~10;                                                             ###
 ### j: AME index 11~12;                                                                         ###
 ### k: data type 0: S2n, 1: S2p; 2: S1n                                                         ###
 ### l: 0: neutron sequence, 1: residue sequence, this index is necessary for python plotting    ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #Dictionary of arrays for plotting purpose, each array key has 5 elements:
 #(proton #; theory masstable index 1~10; AME index 11~12; data type 0: S2n, 1: S2p, 2: S1n; 0: neutron array, 1: residue array)
 ResC = {}
 #Data storage for S2n vs Neutron plot
 # model 1~10, 14,15,17; exp 11~13
 for i in it.chain(range(1,11), range(14,17)):
  for j in range(11,14):
   #FOR NOW WE ONLY PLOT EVEN-Z DATA
    for Z in range(2,zMax+1,2):
     ResC[(Z,i,j,0,0)] = [] #4th index =0: S2n data, 5th index =0: neutron sequence to plot as x-axis
     ResC[(Z,i,j,0,1)] = [] #4th index =0: S2n data, 5th index =1: residue data sequence to plot as y-axis
     #FOR NOW WE ONLY PLOT EVEN-N DATA
     for N in range(2,nMax+1,2):
      if ( (N,Z,i,j) in ResS2n ):
       ResC[(Z,i,j,0,0)].extend([N])
       ResC[(Z,i,j,0,1)].extend([ResS2n[(N,Z,i,j)]])
 #(neutron #; theory masstable index 1~10; AME index 11~12; data type 0: S2n, 1: S2p; 0: neutron sequence, 1: residue sequence)
 #Data storage for S2p vs Proton plot
 # model 1~10, 14,15,16; exp 11~13
 for i in it.chain(range(1,11), range(14,17)):
  for j in range(11,14):
   #FOR NOW WE ONLY PLOT EVEN-N DATA
    for N in range(2,nMax+1,2):
     ResC[(N,i,j,1,0)] = [] #4th index =1: S2p data, 5th index =0: proton sequence to plot as x-axis
     ResC[(N,i,j,1,1)] = [] #4th index =1: S2p data, 5th index =1: residue data sequence to plot as y-axis
     #FOR NOW WE ONLY PLOT EVEN-Z DATA
     for Z in range(2,zMax+1,2):
      if ( (N,Z,i,j) in ResS2p ):
       ResC[(N,i,j,1,0)].extend([Z])
       ResC[(N,i,j,1,1)].extend([ResS2p[(N,Z,i,j)]])

 ###################################################################################################
 ###                                     4.0. PLOT OPTIONS                                       ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 fFormat = "pdf"          #figure format type when saving
 pMod = "global"          #plot mode "test" for testing, "global" or "local" for actual storing
 #pSave, save plot for 0:none; 1:S2n; 2:S2p; 3:S2n & S2p; 4:new nuclei; 5:contour; 6:all
 pSave = 0
 #zMax,nMax is largest proton/neutron # in provided experimental value
 #Range of Z to plot (for S2n)
 zS = 20; zE = zMax           #do not exceed range of zS=2 to zE=zMax
 #Range of N to plot (for S2p)
 nS = 20; nE = nMax          #do not exceed range of nS=2 to nE=nMax
 #Contour plotting, gaussian smoothing or none
 gs_smooth = True
 # gs_range is how far neighboring even-even neutrons and protons to smooth with
 gs_range = 4
 # gs_sigma defaults to 2, defines how spread the gaussian smoothing should be
 gs_sigma = 4; gs_lambda = 0.5 / (gs_sigma**2)
 # Contour color cutoff:
 c_min = -1.5; c_max = 1.5
 ###################################################################################################
 ###                                     4.1. PLOTTING STARTS                                    ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #Calculated chi-square between two sets of data
 s2n_rms = {}
 for i in it.chain(range(1,11), range(14,17)):
  for j in range(11,14):
   for k in range(0,2):
    plt.clf()
    #String for figure name and directory when saved
    fN = "plots/"+pMod+"/" + mtN[j]+"/"
 ###################################################################################################
 ###                                      4.a. S2N PLOTTING                                      ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
    #k = 1, plotting S2p data:
    if ( k==0 and (pSave == 1 or pSave == 3 or pSave == 6) ):
     s2n_rms[(i,j,k)] = 0
     plt.xlabel("Neutron Number",fontsize=14)
     fN = fN + "S2n/" + mtN[i] + "_" + mtN[j] + "_S2n" + "." + fFormat
     count = 0      #counts how many nuclei data point is in the plot
     for Z in range(zS,zE+1,2):
      plt.plot(np.array(ResC[(Z,i,j,k,0)]), np.array(ResC[Z,i,j,k,1]),'r-')
      count = count + len(np.array(ResC[Z,i,j,k,0]))
      for a in range(0,len(np.array(ResC[Z,i,j,k,1]))):
       s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)] + (ResC[Z,i,j,k,1][a])**2
     #Normalizing chi-square with nuclei count
     s2n_rms[(i,j,k)] = math.sqrt(s2n_rms[(i,j,k)]/count)
     #set ticks and range with xx, yy; set y-axis label with ylabel
     if (i == 11 and j == 12):
      xx = np.arange(0,161,20); yy = np.arange(-2,2.1,0.5)
      plt.title("$S_{2n}$:  "+mtN[i]+" vs "+ mtN[j])
      plt.ylabel("$S_{2n,2003}$ - $S_{2n,2016}$ (MeV)",fontsize=14)
      #On the top right corner, print data points and chi-square
      plt.text(115,1.8,"nuclei count: "+str(count) )
      plt.text(115,1.6,"rms = " + str(round( s2n_rms[(i,j,k)],6 )) )
     elif (i == 12 and j == 11):
      xx = np.arange(0,161,20); yy = np.arange(-2,2.1,0.5)
      plt.title("$S_{2n}$:  "+mtN[i]+" vs "+ mtN[j])
      plt.ylabel("$S_{2n,2016}$ - $S_{2n,2003}$ (MeV)",fontsize=14)
      #On the top right corner, print data points and chi-square
      plt.text(115,1.8,"nuclei count: "+str(count) )
      plt.text(115,1.6,"rms = " + str(round( s2n_rms[(i,j,k)],6 )) )
     else:
      xx = np.arange(0,161,20); yy = np.arange(-10,11,2)
      plt.title("$S_{2n}$:  "+mtN[i]+" vs "+ mtN[j])
      plt.ylabel("$S_{2n,th}$ - $S_{2n,exp}$ (MeV)",fontsize=14)
      #On the top right corner, print data points and chi-square
      plt.text(115,8.5,"nuclei count: "+str(count) )
      plt.text(115,7.5,"rms = " + str(round( s2n_rms[(i,j,k)],6 )) )
     plt.xticks(xx); plt.yticks(yy)
     #draw horizontal dashed line of 0 KeV
     plt.axhline(y=0,linestyle='--',color='k')
     #save plot as fN
     plt.savefig(fN,format=fFormat)
     plt.close()
 ###################################################################################################
 ###                                      4.b. S2P PLOTTING                                      ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
    #k = 1, plotting S2p data:
    if ( k==1 and (pSave == 2 or pSave == 3 or pSave == 6) ):
     s2n_rms[(i,j,k)] = 0
     plt.xlabel("Proton Number",fontsize=14)
     fN = fN + "S2p/" + mtN[i] + "_" + mtN[j] + "_S2p" + "." + fFormat
     count = 0      #counts how many nuclei data point is in the plot
     for N in range(nS,nE+1,2):
      plt.plot(np.array(ResC[(N,i,j,k,0)]), np.array(ResC[N,i,j,k,1]),'r-')
      count = count + len(np.array(ResC[N,i,j,k,0]))
      for a in range(0,len(np.array(ResC[Z,i,j,k,1]))):
       s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)] + (ResC[Z,i,j,k,1][a])**2
     #Normalizing chi-square with nuclei count
     s2n_rms[(i,j,k)] = math.sqrt(s2n_rms[(i,j,k)]/count)
     #set ticks and range with xx, yy; set y-axis label with ylabel
     if (i == 11 and j == 12):
      xx = np.arange(0,121,20); yy = np.arange(-8,9,1)
      plt.ylabel("$S_{2p,2003}$ - $S_{2p,2016}$ (MeV)",fontsize=14)
      plt.title("$S_{2p}$:  "+mtN[i]+" vs "+ mtN[j])
      #On the top right corner, print data points and chi-square
      plt.text(90,7.3,"nuclei count: "+str(count) )
      plt.text(90,6.6,"rms = " + str(round( s2n_rms[(i,j,k)],6 )) )
     elif (i == 12 and j == 11):
      xx = np.arange(0,121,20); yy = np.arange(-8,9,1)
      plt.ylabel("$S_{2p,2016}$ - $S_{2p,2003}$ (MeV)",fontsize=14)
      plt.title("$S_{2p}$:  "+mtN[i]+" vs "+ mtN[j])
      #On the top right corner, print data points and chi-square
      plt.text(90,7.3,"nuclei count: "+str(count) )
      plt.text(90,6.6,"rms= " + str(round( s2n_rms[(i,j,k)],6 )) )
     else:
      xx = np.arange(0,121,20); yy = np.arange(-10,11,2)
      plt.ylabel("$S_{2p,th}$ - $S_{2p,exp}$ (MeV)",fontsize=14)
      plt.title("$S_{2p}$:  "+mtN[i]+" vs "+ mtN[j])
      #On the top right corner, print data points and chi-square
      plt.text(90,8.5,"nuclei count: "+str(count) )
      plt.text(90,7.5,"rms = " + str(round( s2n_rms[(i,j,k)],6 )) )
     plt.xticks(xx); plt.yticks(yy)
     #draw horizontal dashed line of 0 KeV
     plt.axhline(y=0,linestyle='--',color='k')
     #save plot as fN
     plt.savefig(fN, format=fFormat)
     plt.close()
 ###################################################################################################
 ###                                 4.c. NEW NUCLEI POINTS PLOT                                 ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################

 newNS2n = []; newZS2n = []; newNS2p = []; newZS2p = []; cNewS2n = 0; cNewS2p = 0;
 newerNS2n = []; newerZS2n = []; newerNS2p = []; newerZS2p = [];
 oldNS2n = []; oldZS2n = []; oldNS2p = []; oldZS2p = []
 for Z in range(2,zMax+1,2):
  for N in range(2,nMax+1,2):
   for k in range(0,2):
    if ( k == 0 and ( (N,Z,k) in newNucS2n) ):
     if ( newNucS2n[(N,Z,k)] == 1 ):
      newNS2n.append(N)
      newZS2n.append(Z)
      cNewS2n = cNewS2n + 1
     elif ( newNucS2n[(N,Z,k)] == -1 ):
      oldNS2n.append(N)
      oldZS2n.append(Z)
     elif ( newNucS2n[(N,Z,k)] == 2 ):
      newerNS2n.append(N)
      newerZS2n.append(Z)
    if ( k == 1 and ( (N,Z,k) in newNucS2n) ):
     if ( newNucS2n[(N,Z,k)] == 1 ):
      newNS2p.append(N)
      newZS2p.append(Z)
      cNewS2p = cNewS2p + 1
     elif ( newNucS2n[(N,Z,k)] == -1 ):
      oldNS2p.append(N)
      oldZS2p.append(Z)
     elif ( newNucS2n[(N,Z,k)] == 2 ):
      newerNS2p.append(N)
      newerZS2p.append(Z)

 #k = 0, plot S2n new nuclei; k = 1 plot S2p new nuclei; k = 2, plot both new nuclei on one plot
 for k in range(0,3):
  plt.clf()
  plt.xlim([0,162]); plt.ylim([0,122])
  xx = np.arange(0,161,20); yy = np.arange(0,121,20)
  plt.xticks(xx); plt.yticks(yy); plt.grid(ls='--', lw = 0.5, alpha = 0.5)
  plt.xlabel("Neutron Number",fontsize=14)
  plt.ylabel("Proton Number",fontsize=14)
  if ( k == 0 and (pSave == 4 or pSave == 6) ):
   plt.plot( np.array(oldNS2n), np.array(oldZS2n), 'o', c='#33FBFF', ms = 3, label = "AME2003" )
   plt.plot( np.array(newNS2n), np.array(newZS2n), 'o', c='#087DF7', ms = 3, label = "AME2016")
   plt.plot( np.array(newerNS2n), np.array(newerZS2n), '*', c='#F70819', ms = 4, label = "2017" )
   fN = "plots/"+pMod+"/newS2n_nuclei.pdf"
   plt.title(r"Experimental S$_{2n}$ data")
   plt.legend(loc = 4)
   plt.savefig(fN, format=fFormat)
   plt.close()
  elif ( k == 1 and (pSave == 4 or pSave == 6) ):
   plt.plot( np.array(newNS2p), np.array(newZS2p), 'o', c='#C94631', ms = 5 )
   plt.title("AME2016 vs AME 2003 new S2p nuclei")
   fN = "plots/"+pMod+"/newS2p_nuclei."+fFormat
   plt.savefig(fN, format=fFormat)
   plt.close()
  #Plot both new S2n/S2p nuclei on one chart
  elif ( k == 2 and (pSave == 4 or pSave == 6) ):
   plt.plot( np.array(newNS2n), np.array(newZS2n), 'o', c='#319FC9', ms = 5, label='S2n' )
   plt.plot( np.array(newNS2p), np.array(newZS2p), 'o', c='#C94631', ms = 5, label='S2p' )
   plt.title("AME2016 vs AME 2003 new S2n/S2p nuclei")
   fN = "plots/"+pMod+"/newS2_nuclei."+fFormat
   plt.legend(numpoints = 1,loc = 4, fontsize = 10)
   plt.savefig(fN, format=fFormat)
   plt.close()

 ###################################################################################################
 ###                               4.d.0 2D CONTOUR PLOTTING DATA                                ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 
 # ContourD: dictionary for contour plotting
 # key is: (theory masstable index 'i'; exp. masstable index 'j'; datatype k; array content: l)
 # For l=0: array of neutron numbers, later to be used on x-axis; l=1: array of proton numbers, later to be used on y-axis, l=3: array of data, for now it's either S2n or S2p residual, depending on k, for pixel coloring
 contourD = {}
 gsN = 2; gsZ = 2;
 # model 1~10, 14,15,17; exp 11~13
 for i in it.chain(range(1,11), range(14,17)):
  for j in range(11,14):
   for k in range(0,1):
    #for each i,j, initialize a new array for storing neutron/proton/data sequence for scipy contour plotting
    contourD[(i,j,k,0)] = []
    contourD[(i,j,k,1)] = []
    contourD[(i,j,k,2)] = []
    for Z in range(2,zMax+1,2):
     for N in range(2,nMax+1,2):
      if not gs_smooth:
       if k == 0:
        if ( (N,Z,i,j) in ResS2n ):
         contourD[(i,j,k,0)].extend([N])
         contourD[(i,j,k,1)].extend([Z])
         contourD[(i,j,k,2)].extend([ResS2n[(N,Z,i,j)]])
       if k == 1:
        if ( (N,Z,i,j) in ResS2p ):
         contourD[(i,j,k,0)].extend([N])
         contourD[(i,j,k,1)].extend([Z])
         contourD[(i,j,k,2)].extend([ResS2p[(N,Z,i,j)]])
    # Gaussian smoothing starts here:
      if gs_smooth:
       # Initialize / reset gaussian summables
       gs_weight_sum = 0
       gs_value = 0
       if k == 0:
        if ( (N,Z,i,j) in ResS2n ):
#         if j == 11: print (f"{Z}, {N}***************")
          # Range of neighboring nuclei
         for p in range(-gs_range, gs_range+1,2):
          gsZ = int(Z + p + 0.000001)
          for q in range(-gs_range, gs_range+1,2):
           gsN = int(N + q + 0.000001)
#           if j == 11:  print (f"{gsZ}, {gsN}")
           if ( (gsN,gsZ,i,j) in ResS2n ):
            # Gaussian average only if neighbor exists, do not apply reflection symmetry
            gs_weight = math.exp( -gs_lambda*(p**2 + q**2) )
            gs_value = gs_value + gs_weight * ResS2n[(gsN,gsZ,i,j)]
            gs_weight_sum = gs_weight_sum + gs_weight
         gs_value = round(gs_value / gs_weight_sum, 6)
         contourD[(i,j,k,0)].extend([N])
         contourD[(i,j,k,1)].extend([Z])
         contourD[(i,j,k,2)].extend([gs_value])
#         if (i == 1 and j == 11):
#          print (str(Z)+ ", " + str(N) + " : "+ str(gs_value))
       if k == 1:
        if ( (N,Z,i,j) in ResS2p ):
          # Range of neighboring nuclei
         for p in range(-gs_range, gs_range+1,2):
          gsZ = int(Z + p + 0.000001)
          for q in range(-gs_range, gs_range+1,2):
           gsN = int(N + q + 0.000001)
           if ( (gsN,gsZ,i,j) in ResS2p ):
            # Gaussian average only if neighbor exists, do not apply reflection symmetry
            gs_weight = math.exp( -gs_lambda*(p**2 + q**2) )
            gs_value = gs_value + gs_weight * ResS2p[(gsN,gsZ,i,j)]
            gs_weight_sum = gs_weight_sum + gs_weight
         gs_value = round(gs_value / gs_weight_sum, 6)
         contourD[(i,j,k,0)].extend([N])
         contourD[(i,j,k,1)].extend([Z])
         contourD[(i,j,k,2)].extend([gs_value])

 ###################################################################################################
 ###                              4.d.1 2D CONTOUR PLOTTING STARTS                               ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 plt.clf()
 # model 1~10, 14,15; exp 11~13
 # Models subplot UNEDF1, DDME2, FRDM2012, SLY4, DDPC1, HFB24, left to right, up to down
 fig, axes = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (30,10))
 gs = gridspec.GridSpec(2, 3, width_ratios=[1]*3, wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.17, right=0.845)
 indLs = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
 ax_ind = (0,0)
 count_axes = -1
 #fN = "plots/"+pMod+"/" + mtN[j] + "/"
 for j in range(11,13):
  count_axes = -1
  for i in it.chain([6,7,14,3,9,15]):
   count_axes = count_axes + 1
   ax_ind = indLs[count_axes]
   for k in range(0,1):
    axes[ax_ind] = plt.subplot(gs[ax_ind])
    #plt.figure(figsize=(10,6),dpi=2400)
    axes[ax_ind].set_xlim([0,162]); axes[ax_ind].set_ylim([0,122])
    xx = np.arange(20,161,20); yy = np.arange(20,121,20)
    axes[ax_ind].set_xticks(xx); axes[ax_ind].set_yticks(yy)
    if ( k==0 and (pSave == 5 or pSave == 6) ):
     # Use dictionary contourD[(i,j,k,l)]
     # Generate data:
     x = np.array(contourD[(i,j,k,0)])
     y = np.array(contourD[(i,j,k,1)])
     z = np.array(contourD[(i,j,k,2)])
#     count = len(z)
#     for a in range(0,len(z)):
#      s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)] + z[a]**2
#     s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)]/count
     axes[ax_ind].scatter(x, y, c=z, marker="s",linewidths=0,s=35
      ,norm = colors.Normalize(vmin=c_min,vmax=c_max,clip=True)
      ,cmap=plt.cm.get_cmap('bwr'))
     if i in [6,7,14]:
      axes[ax_ind].tick_params(axis="x",direction="inout",length = 10)
      axes[ax_ind].set_xticklabels(labels="",color="w")
      #axes[ax_ind].set_xlabel("Neutron Number",fontsize=14)
     if i in [7,14,9,15]:
      axes[ax_ind].tick_params(axis="y",direction="inout",length = 10)
      axes[ax_ind].set_yticklabels(labels="",color="w")
     axes[ax_ind].grid(ls=':',linewidth=0.5,c='#CCCCCC',alpha=0.5)
     axes[ax_ind].axvline(x=20,ls='-',linewidth=1,c='#999999',alpha=0.5)
     axes[ax_ind].text(20, 5, '20', rotation=90)
     axes[ax_ind].axvline(x=28,ls='-',linewidth=1,c='#999999',alpha=0.5)
     axes[ax_ind].text(28, 5, '28', rotation=90)
     axes[ax_ind].axvline(x=50,ls='-',linewidth=1,c='#999999',alpha=0.5)
     axes[ax_ind].text(50, 5, '50', rotation=90)
     axes[ax_ind].axvline(x=82,ls='-',linewidth=1,c='#999999',alpha=0.5)
     axes[ax_ind].text(82, 5, '82', rotation=90)
     axes[ax_ind].axvline(x=126,ls='-',linewidth=1,c='#999999',alpha=0.5)
     axes[ax_ind].text(126, 7, '126', rotation=90)
     axes[ax_ind].set_title( mtN[i],pad=-20.0, fontsize = 15, fontname = "Times New Roman Bold" )
 
  axes[(1,1)].set_xlabel("Neutron Number", fontsize = 18)
  plt.text(-345,140,"Proton Number", rotation = "vertical", fontsize = 18)
  cmap = mpl.cm.get_cmap('bwr')
  ax1 = fig.add_axes([0.87, 0.10, 0.015, 0.8])
  norm = mpl.colors.Normalize(vmin=c_min, vmax=c_max)
  cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
  cb.set_label(r'$S_{2n,exp}$ - $S_{2n,theory}$ (MeV)', fontsize = 18)


  
  if gs_smooth:
   plt.suptitle("S2n Residuals vs "+mtN[j]+", Gaussian smoothing, $\sigma$=" + str(gs_sigma) + ", range="+ str(gs_range), fontsize = 20 )
   plt.savefig("plots/"+mtN[j]+"_6panels_smooth_s_" + str(gs_sigma) + "_r_"+ str(gs_range) +".pdf", format = "pdf")
  elif not gs_smooth:
   plt.suptitle("S2n Residuals vs "+mtN[j] + " raw data", fontsize = 20 )
   plt.savefig("plots/"+mtN[j]+"_6panels_raw_data.pdf", format = "pdf")
  
    
    
#     if not gs_smooth:
#      plt.title("$S_{2n}$:  "+mtN[i]+" vs "+ mtN[j])
#      fN = fN + "S2n/contour/" + mtN[i] + "_" + mtN[j] + "_S2n" + "." + fFormat
#     elif gs_smooth:
#      plt.title("$S_{2n}$:  "+mtN[i]+" vs "+ mtN[j] + " with Gaussian smoothing $\sigma$= "+str(gs_sigma) + ", range= " +str(gs_range) )
#      fN = fN + "S2n/contour/" + mtN[i] + "_" + mtN[j] + "_S2n" + "sigma_" + str(gs_sigma) + ", range= " +str(gs_range) + "." + fFormat
#     plt.text(10,110,"nuclei count: "+str(count) )
   #  plt.text(10,100,"$\chi^2$= " + str(round( s2n_rms[(i,j,k)],6 )) )
#     plt.grid(ls=':',linewidth=0.5,c='#CCCCCC',alpha=0.5)
     #neutron magic number axes
#     plt.axvline(x=20,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(20, 5, '20', rotation=90)
#     plt.axvline(x=28,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(28, 5, '28', rotation=90)
#     plt.axvline(x=50,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(50, 5, '50', rotation=90)
#     plt.axvline(x=82,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(82, 5, '82', rotation=90)
#     plt.axvline(x=126,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(126, 7, '126', rotation=90)
#     plt.savefig(fN, format=fFormat)
#     plt.close()
#    if ( k==1 and (pSave == 5 or pSave == 6) ):
#     # Use dictionary contourD[(i,j,k,l)]
#     # Generate data:
#     x = np.array(contourD[(i,j,k,0)])
#     y = np.array(contourD[(i,j,k,1)])
#     z = np.array(contourD[(i,j,k,2)])
#     count = len(z)
##     for a in range(0,len(z)):
##      s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)] + z[a]**2
##     s2n_rms[(i,j,k)] = s2n_rms[(i,j,k)]/count
#     plt.scatter(x, y, c=z, marker="s",linewidths=0,s=35
#      ,norm = colors.Normalize(vmin=c_min,vmax=c_max,clip=True)
#      ,cmap=plt.cm.get_cmap('bwr'))
#     cb = plt.colorbar()
#     cb.set_label(r'$S_{2p,th}$ - $S_{2p,exp}$ (MeV)')
#     plt.xlabel("Neutron Number",fontsize=14)
#     plt.ylabel("Proton Number",fontsize=14)
#     if not gs_smooth:
#      plt.title("$S_{2p}$:  "+mtN[i]+" vs "+ mtN[j])
#      fN = fN + "S2p/contour/" + mtN[i] + "_" + mtN[j] + "_S2p" + "." + fFormat
#     elif gs_smooth:
#      plt.title("$S_{2p}$:  "+mtN[i]+" vs "+ mtN[j] + " with Gaussian smoothing $\sigma$ = "+str(gs_sigma) + ", range= " +str(gs_range) )
#      fN = fN + "S2p/contour/" + mtN[i] + "_" + mtN[j] + "_S2p" + "sigma_" + str(gs_sigma) + ", range= " +str(gs_range) + "." + fFormat
#     plt.text(10,110,"nuclei count: "+str(count) )
#     #plt.text(10,100,"$\chi^2$= " + str(round( s2n_rms[(i,j,k)],6 )) )
#     plt.grid(ls=':',linewidth=0.5,c='#CCCCCC',alpha=0.5)
#     #proton magic number axes
#     plt.axhline(y=20,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(5, 20, '20')
#     plt.axhline(y=28,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(5, 28, '28')
#     plt.axhline(y=50,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(5, 50, '50')
#     plt.axhline(y=82,ls='-',linewidth=1,c='#999999',alpha=0.5)
#     plt.text(5, 82, '82')
#     #plt.savefig(fN, format=fFormat)
#     plt.close()

 ###################################################################################################
 ###                             5.a. NUMERIC DATA OUTPUT, residual                              ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #Save file option, 0 for none, 1 for residual only, 2 for S2n/p only, 3 for both, 4 for S1n only
 # Do 2:0,2 ; 4:1,3
 tFormat = "csv"
 saveNum = 4; odevity = 1
 
 if ( saveNum == 1 or saveNum == 3):
  outputLabel = "Z".ljust(5) + "N".ljust(5) +  "SKM*(MeV)".ljust(CC) +  "SKP(MeV)".ljust(CC) +  "SLY4(MeV)".ljust(CC) +  "SVMIN(MeV)".ljust(CC) +  "UNEDF0(MeV)".ljust(CC) +  "UNEDF1(MeV)".ljust(CC) +  "UNEDF2(MeV)".ljust(CC) +  "DDME2(MeV)".ljust(CC) +  "DDMEd(MeV)".ljust(CC) +  "DDPC1(MeV)".ljust(CC) +  "NL3*(MeV)".ljust(CC) + "FDRM2012(MeV)".ljust(CC) +  "HFB24(MeV)".ljust(CC)
  outputStr0 = ""; outputStr1 = ""
  #count0, count1:scalars, count how many theories have residue for a specific nuclei
  #if all 10 theories don't have residue, count=0, then don't write that row of Z,N
  count0 = 0; count1 = 0
  #j=11:AME2003; j=12:AME2016
  j = 12
  #remember to change name of output file according to AME data year
  output0 = open("data/All_ResS2n_2016."+tFormat, "w")
  output1 = open("data/All_ResS2p_2016."+tFormat, "w")
  output0.write(outputLabel); output1.write(outputLabel)
  for Z in range(2,zMax,2):
   for N in range(2,nMax,2):
    outputStr0 = "\n" + str(Z).ljust(5) + str(N).ljust(5)
    outputStr1 = "\n" + str(Z).ljust(5) + str(N).ljust(5)
    count0 = 0; count1 = 0
    for i in it.chain(range(1,7), [16], range(7,11), range(14,16)):
     #print (mtN[i])
     #S2n data writing
     if ( (N,Z,i,j) in ResS2n  ):
      outputStr0 = outputStr0 + str(round(ResS2n[(N,Z,i,j)]+0.00000001,6) ).ljust(CC)
      count0 = count0 + 1
     else:
      outputStr0 = outputStr0 + "*".ljust(CC)
     #S2p data writing
     if ((N,Z,i,j) in ResS2p ):
      outputStr1 = outputStr1 + str(round(ResS2p[(N,Z,i,j)]+0.00000001,6) ).ljust(CC)
      count1 = count1 + 1
     else:
      outputStr1 = outputStr1 + "*".ljust(CC)
    #only write to output if at least 1 theory has residual data for (N,Z)
    if (count0 != 0 ):
     output0.write(outputStr0)
    if (count1 != 0 ):
     output1.write(outputStr1)

  output0.close(); output1.close()

 ###################################################################################################
 ###                                  5.b. NUMERIC DATA OUTPUT                                   ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!! Masstable index:                                                                            !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 14-FRDM2012, 15-HFB24, 16-UNEDF2, 17-TRIUMF/RIKEN         !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #S2n
 if ( saveNum == 2 or saveNum == 3):
  CC = 18 #define padding for each colomn
  #outputLabel = "Z".ljust(5) + "N".ljust(5) +  "AME2003(MeV)".ljust(CC) +  "AME2016(MeV)".ljust(CC) +  "JYFLTRAP(MeV)".ljust(CC+1) +  "err2003(MeV)".ljust(CC) +  "err2016(MeV)".ljust(CC) +  "errJYFLTRAP(MeV)".ljust(CC+1) +  "SKM*(MeV)".ljust(CC) +  "SKP(MeV)".ljust(CC) +  "SLY4(MeV)".ljust(CC) +  "SVMIN(MeV)".ljust(CC) +  "UNEDF0(MeV)".ljust(CC) +  "UNEDF1(MeV)".ljust(CC) +  "UNEDF2(MeV)".ljust(CC) + "FRDM2012(MeV)".ljust(CC) +  "HFB24(MeV)".ljust(CC)
  outputLabel = "Z," + "N," +  "2003," +  "2003_sd," +  "2016," +  "2016_sd," +  "2017," +  "2017_sd,"  +  "2018," +  "2018_sd," +\
                "SKM*," +  "SKP," +  "SLY4," +  "SVMIN," +  "UNEDF0," +  "UNEDF1," +  "UNEDF2," + "FRDM2012," +  "HFB24"
  #"DDME2(MeV)".ljust(CC) +  "DDMEd(MeV)".ljust(CC) +  "DDPC1(MeV)".ljust(CC) +  "NL3*(MeV)".ljust(CC)
  outputStr = ""
  #count0, count1:scalars, count how many theories have residue for a specific nuclei
  #if all 10 theories don't have residue, count=0, then don't write that row of Z,N
  count = 0
  #remember to change name of output file according to AME data year
  # even-even S2n only
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/All_S2n_2018_evenZ_evenN."+tFormat, "w")
  # even Z, odd N S2n only
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/All_S2n_2018_evenZ_oddN."+tFormat, "w")
  # odd Z, even N S2n only
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/All_S2n_2018_oddZ_evenN."+tFormat, "w")
  # odd Z, odd N S2n only
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/All_S2n_2018_oddZ_oddN."+tFormat, "w")
  # All S2n available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/All_S2n_2018."+tFormat, "w")
  output.write(outputLabel)
  points = {}
  #print (len(mtN))
  for i in range(1,18):
    points[(mtN[i])] = 0
  for Z in range(zS,zMax1,zOdd):
   for N in range(nS,nMax1,nOdd):
    #outputStr = "\n" + str(Z).ljust(5) + str(N).ljust(5)
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0;
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF/RIKEN 2018
    for i in it.chain(range(11,14), [17]):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S2nErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      outputStr = outputStr + "*" + "," + "*" + ","
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      outputStr = outputStr + "*" + "," #.ljust(CC)
    for i in range(15,16):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      outputStr = outputStr + "*" #.ljust(CC)
    #only write to output if at least 1 set has S2n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr)
  output.close()
 #print (points)
 #S1n
 if ( saveNum == 4):
  #outputLabel = "Z".ljust(5) + "N".ljust(5) +  "AME2003(MeV)".ljust(CC) +  "AME2016(MeV)".ljust(CC) +  "JYFLTRAP(MeV)".ljust(CC+1) +  "err2003(MeV)".ljust(CC) +  "err2016(MeV)".ljust(CC) +  "errJYFLTRAP(MeV)".ljust(CC+1) +  "SKM*(MeV)".ljust(CC) +  "SKP(MeV)".ljust(CC) +  "SLY4(MeV)".ljust(CC) +  "SVMIN(MeV)".ljust(CC) +  "UNEDF0(MeV)".ljust(CC) +  "UNEDF1(MeV)".ljust(CC) + "UNEDF2(MeV)".ljust(CC) + "FRDM2012(MeV)".ljust(CC) +  "HFB24(MeV)".ljust(CC)
  outputLabel = "Z," + "N," +  "2003," +  "2003_sd," +  "2016," +  "2016_sd," +  "2017," +  "2017_sd,"  +  "2018," +  "2018_sd," +\
                "SKM*," +  "SKP," +  "SLY4," +  "SVMIN," +  "UNEDF0," +  "UNEDF1," +  "UNEDF2," + "FRDM2012," +  "HFB24"
  outputStr = ""
  #count0, count1:scalars, count how many theories have residue for a specific nuclei
  #if all theories don't have residue, count=0, then don't write that row of Z,N
  count = 0
  # even-even S1n only
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/All_S1n_2018_evenZ_evenN."+tFormat, "w")
  # even Z, odd N S1n only
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/All_S1n_2018_evenZ_oddN."+tFormat, "w")
  # odd Z, even N S1n only
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/All_S1n_2018_oddZ_evenN."+tFormat, "w")
  # odd Z, odd N S1n only
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/All_S1n_2018_oddZ_oddN."+tFormat, "w")
  # All S1n available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/All_S1n_2018."+tFormat, "w")
  # All S1n available
  elif odevity == 5 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 1;
      output = open("data/test_set_S1n_2018_all_evenZ."+tFormat, "w")
  output.write(outputLabel)

  for Z in range(zS,zMax1,zOdd):
   for N in range(nS,nMax1,nOdd):
    #outputStr = "\n" + str(Z).ljust(5) + str(N).ljust(5)
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF/RIKEN 2018
    for i in it.chain(range(11,14), [17]):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S1nErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
     else:
      outputStr = outputStr + "*" + "," + "*" + ","
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
     else:
      outputStr = outputStr + "*" + "," #.ljust(CC)
    for i in range(15,16):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
     else:
      outputStr = outputStr + "*" #.ljust(CC)
    #only write to output if at least 1 set has S1n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr)
  output.close()

 ###################################################################################################
 ###                           6. DRIPLINES & S2N AVAILABLE REGIONS                              ###
 ###################################################################################################
 # Label nuclei with S2n or S2p data in 2016 with 

 ###################################################################################################
 ###                                      FUNCTION EXECUTION                                     ###
 ###################################################################################################



residues()

