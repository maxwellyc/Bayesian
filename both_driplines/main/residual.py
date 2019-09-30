###################################################################################################
### This Code was written to extract and plot residuals of 2 nucleon separation energy from     ###
### theoretical masstable calculations with Skyrme DFT vs experimental data                     ###
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.mlab import bivariate_normal
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm
import scipy.interpolate
import scipy.ndimage as ndimage
import itertools as it
import math


# For Adobe illustrator text
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#For calculation of odd values please refer to: http://massexplorer.frib.msu.edu/content/masstables/Odd_Values_from_Even_Data.pdf

#mtN stands for masstableName
mtN = ["","SkM*","SkP","SLy4","SV-min","UNEDF0","UNEDF1","DD-ME2","DD-ME$\delta$","DDPC1","NL3*", #1~10
        "AME2003","AME2016","JYFLTRAP2017","FRDM2012","HFB24","UNEDF2","TRIUMF2018","RIKEN2018","new_other", #11~19
        "","","","oct-SkM*","oct-SkP","oct-SLy4","oct-SV-min","oct-UNEDF0","oct-UNEDF1","oct-UNEDF2"] #23~29

#low_lim is a filter for separation value, eg. when low_lim is 0 only non-negative ( >= 0) values will be recorded
def data_import(low_lim = -100000):
 # function that detects if arg. is a number

 # Atomic mass unit from Wikipedia
 # Note that FRDM2012 uses 1u = 931.5014MeV
 uAtomic = 931.49        # MeV rounded to 2 decimals
 uAtomic0 = 931494.0954  # keV
 # uAtomic0 = 931494.061 # Erik's uAtomic, keV

 # Hydrogen atom and neutron masses From AME 2012, which HFB-24 was fitted to
 # data before 2/13/19 were computed by these values
 M_P0 = 1.0078250322 * uAtomic0  #938.783 MeV
 M_N0 = 1.0086649158 * uAtomic0  #939.565 MeV
 M_P = 938.78  # MeV rounded to 2 decimals  # used M_P = 938.27 before 2/13/19 incorrectly
 M_N = 939.57  # MeV rounded to 2 decimals
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

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
 dir = "data/max_raw_no_dripline/"; tail = "_all_nuclei_max.dat"
 f1 = open(dir+"SKMS"+tail); f2 = open(dir+"SKP"+tail); f3 = open(dir+"SLY4"+tail)
 #f3 = open(dir+"SLY4_massexplorer_separation.dat"); # testing massexplorer data
 #f3 = open(dir+"SLY4_octupole_separation_filtered.dat") # testing octupole residule data
 f4 = open(dir+"SV-MIN"+tail); f5 = open(dir+"UNEDF0"+tail)
 f6 = open(dir+"UNEDF1"+tail); f16 = open(dir+"UNEDF2"+tail)
 #f16 = open(dir+"UNEDF2_octupole_separation.dat") # testing octupole residule data
 #column indices of Z, N, S2n, S2p, S1n, S1p, binding energy. ME stands for MassExplorer
 aME = [1,2,8,6,7,5,4]

 # 2019 HFBTHO octupole mass-tables:
 # refined_ground_state_selection disregards energy that corresponds to beta2/3 > 0.5
 #
 refined_gs = True
 if refined_gs:
    dir2 = "data/octupole_2019/refined/"; tail2 = "_binding_beta2_3_filtered_refined.dat"
 elif not refined_gs:
    dir2 = "data/octupole_2019/unrefined/"; tail2 = "_binding_beta2_3_filtered.dat"
 f23 = open(dir2+"SKMS"+tail2); f24 = open(dir2+"SKP"+tail2); f25 = open(dir2+"SLY4"+tail2)
 f26 = open(dir2+"SV-MIN"+tail2); f27 = open(dir2+"UNEDF0"+tail2)
 f28 = open(dir2+"UNEDF1"+tail2); f29 = open(dir2+"UNEDF2"+tail2)
 #column indices of Z, N, binding energy, beta2_out, beta3_out, OCT stands for octupole
 aOCT = [0,1,3,4,5]


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

 #2003 & 2016 AME separation energies
 f11 = open("data/2003AME_S2.dat"); f12 = open("data/2016AME_S2.dat")
 f17= open("data/2003AME_S1.dat"); f18 = open("data/2016AME_S1.dat")
 #AME stands for Atomic Mass Evaluation, this aAME is good for both *_S2 file and *_S1 file
 #column indices of Z, N, S2n/S1n, S2n/S1n errors, S2p/S1p, S2p/S1p errors (error data only exists for AME)
 aAME = [0,1,2,3,4,5]

 #2016 AME mass excess raw input
 #Z,N, mass excess, mass excess error
 f22 = open("data/2016AME_ME.dat")
 aAME0 = [1,0,5,6]

 #Additional experimental data after 2016:
 # JYFLTRAP 2017 rare earth see email communications w/ Witek
 # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.262701  158Nd; 160Pm; 162Sm; 164-166Gd; JYFLTRAP2017
 # https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.062503       Ti, TRIUMF 2018
 # https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022506       55-57Ca, RIKEN 2018
 # https://journals.aps.org/prc/abstract/10.1103/PhysRevC.96.014310      101,102Rb
 # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.152501  246Es; 249,250,252Md
 # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.262702  159-160Nd; 163-164Sm. others included in JYFLTRAP2017
 # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.192502  77,79Cu
 f13 = open("data/JYFL2017.dat"); f19 = open("data/TRIUMF2018.dat"); f20 = open("data/RIKEN2018.dat")
 f21 = open("data/new_other.dat")
 #column indices of Z, N, mass excess(keV), mass excess error(keV)
 aJYFL = [0,1,4,5]; a2018 = [0,1,2,3]

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
 # Mass Excess from 2018 TRIUMF / RIKEN / new other / AME2016
 l[19] = f19.readlines(); l[20] = f20.readlines(); l[21] = f21.readlines(); l[22] = f22.readlines()
 # 2019 Octupole tables:
 l[23] = f23.readlines(); l[24] = f24.readlines(); l[25] = f25.readlines(); l[26] = f26.readlines();
 l[27] = f27.readlines(); l[28] = f28.readlines(); l[29] = f29.readlines()

 f1.close(); f2.close(); f3.close(); f4.close(); f5.close(); f6.close(); f7.close(); f8.close()
 f9.close(); f10.close(); f11.close(); f12.close(); f13.close(); f14.close(); f15.close()
 f16.close(); f17.close(); f18.close(); f19.close(); f20.close(); f21.close(); f22.close()
 f23.close(); f24.close(); f25.close(); f26.close(); f27.close(); f28.close(); f29.close()
 ###################################################################################################
 ### S2* dictionaries contain all the separation energies read from the 12 masstables above      ###
 ### The key is structured as: (N,Z,masstable index),                                            ###
 ### Error dictionaries only contain errors from AME2003 and AME2016                             ###
 ### To aviod confusion, the masstable index of error dictionaries for AME2003/AME2016           ###
 ### is still 11/12, in the future it's possible to include errors for theories                  ###
 ### (Neutron #, Proton #, theory masstable index, exp. masstable index )                        ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

 ###################################################################################################
 ###                                  2. READ IN DATA FROM FILE                                  ###
 ###################################################################################################
 S2n = {}; S2nErr = {}; S2p = {}; S2pErr = {}; BE={}; BEErr={}; nMax = 2; zMax = 2 ; nMax1 = 2; zMax1 = 2
 S1n = {}; S1nErr = {}; S1p = {}; S1pErr = {}
 ###################################################################################################
 ###                            2.a Mass Explorer SKYRME ::READ-INPUT                            ###
 ###################################################################################################
 #Mass Explorer data read and store, use aME:
 for i in it.chain(range(1,7), range(16,17)):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(ss[int(aME[1])])      #Neutron Number
     Z = int(ss[int(aME[0])])      #Proton  Number
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
     #Binding Energy
     data = str(ss[int(aME[6])])
     if (isNum(data) and float(data) < 0):  # there's a filter of requiring BE to be negative, o/w pointless
      BE[(N,Z,i)] = float(data)
    except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue

# Find abnormal S2n values
# for i in it.chain(range(1,7), range(16,17)):
#  for Z in range(2,zMax1,2):
#    for N in range(126,nMax1,2):
#        if (N,Z,i) in S2n:
#            if S2n[(N,Z,i)] > 30 :
#                print (mtN[i],Z,N,S2n[(N,Z,i)],S2p[(N,Z,i)])
 ###################################################################################################
 ###                        2.b Relativistic Mean Field ::READ-INPUT                             ###
 ###################################################################################################
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
 ###################################################################################################
 ###                          2.c Atomic Mass Evaluations ::READ-INPUT                           ###
 ###################################################################################################
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
     if (isNum(data) and float(data) >= low_lim ):
      S2n[(N,Z,i)] = float(data)#/1000.0             #If experimental data is in keV, convert to MeV
      #S2n Error, there's an indent here, only store error if there's S2n data
      data = str(ss[int(aAME[3])])
      if (isNum(data)):
       S2nErr[(N,Z,i)] = float(data)#/1000.0         #If experimental data is in keV, convert to MeV
     #S2p
     data = str(ss[int(aAME[4])])
     #print (mtN[i],Z,N,data)
     if (isNum(data) and float(data) >= low_lim ):
      S2p[(N,Z,i)] = float(data)#/1000.0             #If experimental data is in keV, convert to MeV
      #S2p Error, there's an indent here, only store error if there's S2p data
      data = str(ss[int(aAME[5])])
      if (isNum(data)):
       S2pErr[(N,Z,i)] = float(data)#/1000.0         #If experimental data is in keV, convert to MeV
    except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue
 #S1n & S1p:
 for i in range(17,19):
  for line in l[i]:
    # This i to ii ensure that S1n[(N,Z,11)] is from AME2003, S1n[(N,Z,12)] from AME2016, to avoid confusion
    # The difference in index for experimental S1n is contained in the data reading stage only
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
     if (isNum(data) and float(data) >= low_lim ):
      S1n[(N,Z,ii)] = float(data)#/1000.0             #If experimental data is in keV, convert to MeV
      #S1n Error, there's an indent here, only store error if there's S1n data
      data = str(ss[int(aAME[3])])
      if (isNum(data)):
       S1nErr[(N,Z,ii)] = float(data)#/1000.0         #If experimental data is in keV, convert to MeV
     data = str(ss[int(aAME[4])])
     if (isNum(data) and float(data) >= low_lim ):
      S1p[(N,Z,ii)] = float(data)#/1000.0             #If experimental data is in keV, convert to MeV
      #S1n Error, there's an indent here, only store error if there's S1n data
      data = str(ss[int(aAME[5])])
      if (isNum(data)):
       S1pErr[(N,Z,ii)] = float(data)#/1000.0         #If experimental data is in keV, convert to MeV
    except (ValueError, IndexError):        #N,Z, or, S1n are not numbers
       continue
 #AME2016 binding energy from mass excess
 for i in range(12,13):
  for line in l[22]:
    ss = line.split("\t")
    try:
     N = int(float(ss[int(aAME0[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aAME0[0])])+0.0001)      #Proton Number
     A = N+Z
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     #mass excess and error
     data = str(ss[int(aAME0[2])])
     if (isNum(data)):
      BE[(N,Z,i)] = round(-0.001 * ( Z*M_P0 + N*M_N0 - (float(data) + A*uAtomic0) ), 6) # Binding energy set to negative
      data = str(ss[int(aAME0[3])])
      if (isNum(data)):
       BEErr[(N,Z,i)] = 0.001 * float(data)
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 ###################################################################################################
 ###                      2.d New experimental data after 2016 ::READ-INPUT                      ###
 ###################################################################################################
 #JYFLTRAP 2017 binding energy from mass excess
 for i in range(13,14):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aJYFL[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aJYFL[0])])+0.0001)      #Proton Number
     A = N+Z
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     #mass excess and error
     data = str(ss[int(aJYFL[2])])
     if (isNum(data)):
      BE[(N,Z,i)] = round(-0.001 * ( Z*M_P0 + N*M_N0 - (float(data) + A*uAtomic0) ), 6)
      data = str(ss[int(aJYFL[3])])
      if (isNum(data)):
       BEErr[(N,Z,i)] = 0.001 * float(data)
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 #New data (TRIUMF, RIKEN, other newer data), use a2018
 for i in range(19,22):
  for line in l[i]:
    if (i == 19) : ii = 17
    if (i == 20) : ii = 18
    if (i == 21) : ii = 19
    ss = line.split()
    try:
     N = int(float(ss[int(a2018[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(a2018[0])])+0.0001)      #Proton Number
     A = N+Z
     if (nMax<N):
      nMax = N
     if (zMax<Z):
      zMax = Z
     #mass excess and error
     data = str(ss[int(a2018[2])])
     if (isNum(data)):
      BE[(N,Z,ii)] = round(-0.001 * ( Z*M_P0 + N*M_N0 - (float(data) + A*uAtomic0) ), 6)
      data = str(ss[int(a2018[3])])
      if (isNum(data)):
       BEErr[(N,Z,ii)] = 0.001 * float(data)
    except (ValueError, IndexError):                #N,Z, or, S2n/S2p are not numbers
       continue
 ###################################################################################################
 ###                          2.e Phenomenological Models ::READ-INPUT                           ###
 ###################################################################################################

 #FRDM 2012, USE aFRDM
 for i in range(14,15):
  for line in l[i]:
    ss = line.split()
    try:
     N = int(float(ss[int(aFRDM[1])])+0.0001)      #Neutron Number
     Z = int(float(ss[int(aFRDM[0])])+0.0001)      #Proton Number
     #S2n
     if (nMax1<N):
      nMax1 = N
     if (zMax1<Z):
      zMax1 = Z
     data = str(ss[int(aFRDM[2])])
     if (isNum(data)):
      BE[(N,Z,i)] = -1.0 * float(data)              #Binding energy is positive in FRDM2012 file, we use negative convention
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
     if (nMax1<N):
      nMax1 = N
     if (zMax<Z):
      zMax1 = Z
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
 ###################################################################################################
 ###                 2.f  Calculate separation energies from HFB-24 & FRDM2012                   ###
 ###################################################################################################
 # calulate S2n/p for HFB-24 & FRDM2012
 d_cache = 0.0
 for i in range(14,16):
  for Z in range(2,zMax1+1):
   for N in range(2,nMax1+1):
    if ( ( (N,Z,i) in BE) and ( (N-2,Z,i) in BE) ):
     d_cache =  - BE[(N,Z,i)] + BE[(N-2,Z,i)]
     if d_cache >= low_lim : S2n[(N,Z,i)] = d_cache
    if ( ( (N,Z,i) in BE) and ( (N,Z-2,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-2,i)]
     if d_cache >= low_lim : S2p[(N,Z,i)] = d_cache

 # calculate S1n/S1p for FRDM2012, HFB-24's are read-in
 for i in range(14,15):
  for Z in range(2,zMax1+1):
   for N in range(2,nMax1+1):
    if ( ( (N,Z,i) in BE) and ( (N-1,Z,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N-1,Z,i)]
     if  d_cache >= low_lim : S1n[(N,Z,i)] = d_cache
    if ( ( (N,Z,i) in BE) and ( (N,Z-1,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-1,i)]
     if  d_cache >= low_lim : S1p[(N,Z,i)] = d_cache

 ###################################################################################################
 ###                  2.g  Calculate updated experimental separation energies                    ###
 ###################################################################################################
 #!! 11-AME2003, 12-AME2016, 13-JYFL2017, 17-TRIUMF, 18-RIKEN, 19-new_other  !!#
 # If 13, 17-19 has data points, do a search in its vicinity, any S1/2n, S1/2p generated using it should be stored in 13, 17-19
 for i in it.chain([13,17,18,19]):
  for Z in range(2,zMax+1):
   for N in range(2,nMax+1):
    if (N,Z,i) in BE:
     ### NEUTRON SEPARATION ###
     ## S1n of up and down, if mass of neighboring nuclei not in new dataset, check AME2016
     if (N-1,Z,i) in BE:
      S1n[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N-1,Z,i)], 6)
      S1nErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N-1,Z,i)]**2.0), 6)
     elif (N-1,Z,12) in BE and (N-1,Z,i) not in BE: #2016
      S1n[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N-1,Z,12)], 6)
      S1nErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N-1,Z,12)]**2.0), 6)
     if (N+1,Z,i) in BE:
      S1n[(N+1,Z,i)] = -1.0 * round( BE[(N+1,Z,i)] - BE[(N,Z,i)], 6)
      S1nErr[(N+1,Z,i)] = round(math.sqrt(BEErr[(N+1,Z,i)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     elif (N+1,Z,12) in BE and (N+1,Z,i) not in BE: #2016
      S1n[(N+1,Z,i)] = -1.0 * round( BE[(N+1,Z,12)] - BE[(N,Z,i)], 6)
      S1nErr[(N+1,Z,i)] = round(math.sqrt(BEErr[(N+1,Z,12)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     ## S2n of up and down, if mass of neighboring nuclei not in new dataset, check AME2016
     if (N-2,Z,i) in BE:
      S2n[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N-2,Z,i)], 6)
      S2nErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N-2,Z,i)]**2.0), 6)
     elif (N-2,Z,12) in BE and (N-2,Z,i) not in BE: #2016
      S2n[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N-2,Z,12)], 6)
      S2nErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N-2,Z,12)]**2.0), 6)
     if (N+2,Z,i) in BE:
      S2n[(N+2,Z,i)] = -1.0 * round( BE[(N+2,Z,i)] - BE[(N,Z,i)], 6)
      S2nErr[(N+2,Z,i)] = round(math.sqrt(BEErr[(N+2,Z,i)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     elif (N+2,Z,12) in BE and (N+2,Z,i) not in BE: #2016
      S2n[(N+2,Z,i)] = -1.0 * round( BE[(N+2,Z,12)] - BE[(N,Z,i)], 6)
      S2nErr[(N+2,Z,i)] = round(math.sqrt(BEErr[(N+2,Z,12)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)

     ### PROTON SEPARATION ###
     ## S1p of up and down, if mass of neighboring nuclei not in new dataset, check AME2016
     if (N,Z-1,i) in BE:
      S1p[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N,Z-1,i)], 6)
      S1pErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N,Z-1,i)]**2.0), 6)
     elif (N,Z-1,12) in BE and (N,Z-1,i) not in BE: #2016
      S1p[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N,Z-1,12)], 6)
      S1pErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N,Z-1,12)]**2.0), 6)
     if (N,Z+1,i) in BE:
      S1p[(N,Z+1,i)] = -1.0 * round( BE[(N,Z+1,i)] - BE[(N,Z,i)], 6)
      S1pErr[(N,Z+1,i)] = round(math.sqrt(BEErr[(N,Z+1,i)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     elif (N,Z+1,12) in BE and (N,Z+1,i) not in BE: #2016
      S1p[(N,Z+1,i)] = -1.0 * round( BE[(N,Z+1,12)] - BE[(N,Z,i)], 6)
      S1pErr[(N,Z+1,i)] = round(math.sqrt(BEErr[(N,Z+1,12)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     ## S2p of up and down, if mass of neighboring nuclei not in new dataset, check AME2016
     if (N,Z-2,i) in BE:
      S2p[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N,Z-2,i)], 6)
      S2pErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N,Z-2,i)]**2.0), 6)
     elif (N,Z-2,12) in BE and (N,Z-2,i) not in BE: #2016
      S2p[(N,Z,i)] = -1.0 * round( BE[(N,Z,i)] - BE[(N,Z-2,12)], 6)
      S2pErr[(N,Z,i)] = round(math.sqrt(BEErr[(N,Z,i)]**2.0 + BEErr[(N,Z-2,12)]**2.0), 6)
     if (N,Z+2,i) in BE:
      S2p[(N,Z+2,i)] = -1.0 * round( BE[(N,Z+2,i)] - BE[(N,Z,i)], 6)
      S2pErr[(N,Z+2,i)] = round(math.sqrt(BEErr[(N,Z+2,i)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)
     elif (N,Z+2,12) in BE and (N,Z+2,i) not in BE: #2016
      S2p[(N,Z+2,i)] = -1.0 * round( BE[(N,Z+2,12)] - BE[(N,Z,i)], 6)
      S2pErr[(N,Z+2,i)] = round(math.sqrt(BEErr[(N,Z+2,12)]**2.0 + BEErr[(N,Z,i)]**2.0), 6)

  # Remove duplicates of JYFLTRAP2017 within new_other
  for Z in range(2,zMax+1):
   for N in range(2,nMax+1):
    if (N,Z,13) in S1n and (N,Z,19) in S1n:
     if S1n[(N,Z,13)] == S1n[(N,Z,19)]: del S1n[(N,Z,19)]
    if (N,Z,13) in S1p and (N,Z,19) in S1p:
     if S1p[(N,Z,13)] == S1p[(N,Z,19)]: del S1p[(N,Z,19)]
    if (N,Z,13) in S2n and (N,Z,19) in S2n:
     if S2n[(N,Z,13)] == S2n[(N,Z,19)]: del S2n[(N,Z,19)]
    if (N,Z,13) in S2p and (N,Z,19) in S2p:
     if S2p[(N,Z,13)] == S2p[(N,Z,19)]: del S2p[(N,Z,19)]

 ###################################################################################################
 ###                2.h Octupole Masstable 2019 HFBTHO Maxwell Cao ::READ-INPUT                  ###
 ###################################################################################################
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 ###################################################################################################

 #column indices of Z, N, binding energy, beta2_output, beta3_output, OCT stands for octupole
 # aOCT = [0,1,3,4,5]
 # Mass Explorer data read and store, use aME:
 # Deformation dictionary, key is N,Z, octuple masstable index, value is tuple: (beta2_output, beta3_output)
 deform = {}
 # oct_nuc dictionary, key is N,Z, octuple masstable index, value is boolean, only True for nuclei beta3 >= beta3_cut
 # If you don't need a filter when calculating rms error, use beta3_cut = -1, beta2_cap = beta3_cap = 1000
 beta3_cut = -1; beta3_cap = beta2_cap = 1000
 oct_nuc = {}; oct_nuc1 = {}
 for i in range(23,30):
    for line in l[i]:
        ss = line.split()
        try:
            N = int(ss[int(aOCT[1])]); Z = int(ss[int(aOCT[0])])
            if (nMax1<N): nMax1 = N
            if (zMax1<Z): zMax1 = Z
            # Binding energy
            data = str(ss[int(aOCT[2])])
            if isNum(data) and float(data) < 0 :  # there's a filter of requiring BE to be negative, o/w pointless
                BE[(N,Z,i)] = float(data)
            # Deformation
            data1 = str(ss[int(aOCT[3])]) #beta2
            data2 = str(ss[int(aOCT[4])]) #beta3
            if isNum(data1) and isNum(data2) :
                deform[(N,Z,i)] = (float(data1), float(data2))
                # set plausible range for beta3, although not all nuclei has experimental value
                if (abs(float(data2)) >= beta3_cut and abs(float(data2)) <= beta3_cap
                    and abs(float(data1)) <= beta2_cap):
                    oct_nuc[(N,Z,i)] = True
                else:
                    oct_nuc[(N,Z,i)] = False
                #if abs(float(data2))>=0.5 or abs(float(data1))>=0.5:
                    #print (mtN[i],"\t",Z,"\t",N,"\t",float(data1),"\t",float(data2))
        except (ValueError, IndexError):
            continue
 # check deformation in Z = 120 chain
# for Z in range(120,121,2):
#  for N in range(2,301,2):
#    try:
#        print (Z,N,deform[(N,Z,29)])
#    except:
#        continue
 # set oct_nuc of old masstable data to be true so that we can compare later
 for (N,Z,i) in oct_nuc.keys():
    if 22 < i < 29: oct_nuc1[(N,Z,i-22)] = oct_nuc[(N,Z,i)]
    if i == 29: oct_nuc1[(N,Z,16)] = oct_nuc[(N,Z,i)]
 #merge 2 octupole deformed nuclei list
 oct_nuclei = {**oct_nuc, **oct_nuc1}
 ###################################################################################################
 ###                  2.i Octupole Masstable separation energies calculation                     ###
 ###################################################################################################
 # low_lim sets cutoff of separation energies, sometimes we only need positive ones.
 for i in range(23,30):
  for Z in range(2,zMax1+1):
   for N in range(2,nMax1+1):
    d_cache = 0.0
    if ( ( (N,Z,i) in BE) and ( (N-2,Z,i) in BE) ):
     d_cache =  - BE[(N,Z,i)] + BE[(N-2,Z,i)]
     if d_cache >= low_lim : S2n[(N,Z,i)] = round(d_cache,6)
    if ( ( (N,Z,i) in BE) and ( (N,Z-2,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-2,i)]
     if d_cache >= low_lim : S2p[(N,Z,i)] = round(d_cache,6)
    if ( ( (N,Z,i) in BE) and ( (N-1,Z,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N-1,Z,i)]
     if  d_cache >= low_lim : S1n[(N,Z,i)] = round(d_cache,6)
    if ( ( (N,Z,i) in BE) and ( (N,Z-1,i) in BE) ):
     d_cache = - BE[(N,Z,i)] + BE[(N,Z-1,i)]
     if  d_cache >= low_lim : S1p[(N,Z,i)] = round(d_cache,6)
 ###################################################################################################
 ###                             2.j Dripline evaluation for Skyrme                              ###
 ###################################################################################################
 # rough cutoff energy for dripline, for Bayesian purpose sometimes we don't use exact 0 MeV as cutoff
 # we do a rough even-even cut here, check for consecutive separations below cut to actually perform cut
 s2n_cut = -1; s2p_cut = -1; dripline = {}
 zz = 118; nn = 134; iii = 1
# for N in range(2,nMax1+1,2):
#    if (N,zz,iii) in S2n:
#        print (mtN[iii],zz,N, S2n[(N,zz,iii)], S2p[(N,zz,iii)])
# for Z in range(2,zMax1+1,2):
#    if (nn,Z,iii) in S2p:
#        print (mtN[iii],Z,nn, S2p[(nn,Z,iii)])


 # dripline[(Z,i)]
 for i in it.chain(range(1,7), [16], range(23,30)):
    # Neutron dripline
    for Z in range(2,zMax1+1,2):
        for N in range(2, nMax1+1,2):
            if (N,Z,i) in S2n and (N+2,Z,i) in S2n:
                if S2n[(N,Z,i)] >= s2n_cut and S2n[(N+2,Z,i)] <= s2n_cut and N/Z > 1.5:
                    dripline[(Z,i,1)] = N
                    break
    # Proton dripline
    for Z in range(2,zMax1+1,2):
        for N in reversed(range(2,nMax1+1,2)):
            if (N,Z,i) in S2p and (N-2,Z,i) in S2p:
                if S2p[(N,Z,i)] >= s2p_cut and S2p[(N-2,Z,i)] <= s2p_cut and N/Z < 1.5:
                    dripline[(Z,i,0)] = N
                    break

    # Catch dripline if already filtered
    for Z in range(2,zMax1+1,2):
        for N in range(2, nMax1+1,2):
            if (N,Z,i) in BE and (Z,i,0) not in dripline:
                dripline[(Z,i,0)] = N
                break
        for N in reversed(range(2,nMax1+1,2)):
            if (N,Z,i) in BE and (Z,i,1) not in dripline:
                dripline[(Z,i,1)] = N
                break
 return S2n,S2p,S1n,S1p,BE,oct_nuclei,zMax,nMax
# for i in range(23,30):
#  for Z in range(2,121,2):
#    for N in range(2,301,2):
#        if (N,Z,i) in deform:
#            if abs(deform[(N,Z,i)][0]) > 0.5 or abs(deform[(N,Z,i)][1]) > 0.5:
#                print (mtN[i],"\t",Z,"\t",N,"\t",deform[(N,Z,i)][0],"\t",deform[(N,Z,i)][1],"\t",
#                    S2n[(N,Z,i)] if (N,Z,i) in S2n else "None","\t",
#                    S2p[(N,Z,i)] if (N,Z,i) in S2p else "None","\t",
#                    BE[(N,Z,i)] if (N,Z,i) in BE else "None")
#    for Z in range(2,zMax1+1,2):
#        if (Z,i,0) in dripline and i == 1:
#            print (mtN[i],"\t",Z,"\t",dripline[(Z,i,0)],"\t",dripline[(Z,i,1)])





 ###################################################################################################
 ###                    3. CALCULATE residualS, RESULTS STORED IN ResS2n, ResS2p                  ###
 ###################################################################################################

 ###################################################################################################
 ### ResS2n, ResS2p are Dictionaries for S2n, S2p residuals between theory and exp.               ###
 ### The key is structrued as (N,Z,i,j):                                                         ###
 ### (Neutron #, Proton #, theory masstable index, exp. masstable index )                        ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

def residuals():
 S2n,S2p,S1n,S1p,BE,oct_nuclei,zMax,nMax = data_import()
 ResS2n = {}; ResS2p = {}; ResS1n = {}; ResS1p = {}; ResBE = {}
 rms_s2n = {}; rms_s1n = {}; rms_s2p = {}; rms_s1p = {}; rms_be = {}
 zs = 2; ns = 2
 #res1_p = {}; res2_p = {}
 del_switch = 0 # =1 / 2 to delete from del_list1 / del_list2
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

 # model 1~6, 14~16, 23~29 ; exp 12
 for i in it.chain(range(1,7), range(14,17), range(23,30)):
  #res2_p[mtN[i]] = 0
  # Only AME2016 has binding energy in this code, for now let's just focus on AME2016 residuals
  for j in range(12,13):
    for Z in range(2,zMax+1):
     for N in range(2,nMax+1):
      # S2n residual calculation
      if (N,Z,i) in S2n and (N,Z,j) in S2n :
       ResS2n[(N,Z,i,j)] = round(S2n[(N,Z,j)] - S2n[(N,Z,i)],6)
       #if j == 11 and (not N%2) and (not Z%2): res2_p[mtN[i]] +=  1
      # S2p residual calculation
      if (N,Z,i) in S2p and (N,Z,j) in S2p :
       ResS2p[(N,Z,i,j)] = round(S2p[(N,Z,j)] - S2p[(N,Z,i)],6)
      # S1n residual calculation
      if (N,Z,i) in S1n and (N,Z,j) in S1n :
       ResS1n[(N,Z,i,j)] = round(S1n[(N,Z,j)] - S1n[(N,Z,i)],6)
      # S1p residual calculation
      if (N,Z,i) in S1p and (N,Z,j) in S1p :
       ResS1p[(N,Z,i,j)] = round(S1p[(N,Z,j)] - S1p[(N,Z,i)],6)
      # Binding energy residual calculation
      if (N,Z,i) in BE and (N,Z,j) in BE :
       ResBE[(N,Z,i,j)] = round(BE[(N,Z,j)] - BE[(N,Z,i)],6)

 #Only Skyrme and FRDM-2012, HFB-24 has S1n/p values, RMF's are not calculated

 # Find rms of S2n points in AME2016 but not in AME2003, starting from zs, ns
# for i in it.chain(range(1,7), range(14,17)):
#  rms_s2n[(i,12)] = 0
#  rms_s1n[(i,12)] = 0
#  count_rms = 0
#  count_rms1= 0
#  #S2n
#  for Z in range(zs,zMax+1,2):
#   for N in range(ns,nMax+1,2):
#    if (N,Z,i,12) in ResS2n.keys() and (N,Z,i,11) not in ResS2n :
#     rms_s2n[(i,12)] = rms_s2n[(i,12)] + ResS2n[(N,Z,i,12)]**2
#     count_rms += 1
#  if count_rms != 0:
#    rms_s2n[(i,12)] = math.sqrt(rms_s2n[(i,12)] / float(count_rms))
#
#  #S1n
#  for Z in range(zs,zMax+1):
#   for N in range(ns,nMax+1):
#    if (N,Z,i,12) in  ResS1n and (N,Z,i,11) not in ResS1n :
#     rms_s1n[(i,12)] = rms_s1n[(i,12)] + ResS1n[(N,Z,i,12)]**2
#     count_rms1 = count_rms1 + 1
#  if count_rms1 != 0: rms_s1n[(i,12)] = math.sqrt(rms_s1n[(i,12)] / float(count_rms1))

 # theory rms error when compared with AME2016 with Z>=zs, N>=zs
 # check even,even first, Z, N range has step size 2
 #S2n
 for i in it.chain(range(1,7),[16], range(23,30)):
  rms_s2n[(i,12)] = 0; rms_s2p[(i,12)] = 0; rms_s1n[(i,12)] = 0; rms_s1p[(i,12)] = 0; rms_be[(i,12)] = 0;
  count_rms1 = count_rms2 = count_rms3 = count_rms4 = count_rms5 = 0
  for Z in range(zs,zMax+1,2):
   for N in range(ns,nMax+1,2):
    if (N,Z,i) in oct_nuclei and oct_nuclei[(N,Z,i)]: # Use this statement to compare only octupole deformed nuclei
        if (N,Z,i,12) in ResS2n:
            rms_s2n[(i,12)] += ResS2n[(N,Z,i,12)]**2
            count_rms1 += 1
        if (N,Z,i,12) in ResS2p and not (Z == 4 and N == 8):
            rms_s2p[(i,12)] += ResS2p[(N,Z,i,12)]**2
            count_rms2 += 1
        if (N,Z,i,12) in ResS1n:
            rms_s1n[(i,12)] += ResS1n[(N,Z,i,12)]**2
            count_rms3 += 1
        if (N,Z,i,12) in ResS1p:
            rms_s1p[(i,12)] += ResS1p[(N,Z,i,12)]**2
            count_rms4 += 1
        if (N,Z,i,12) in ResBE:
            rms_be[(i,12)] += ResBE[(N,Z,i,12)]**2
            count_rms5 += 1
  rms_s2n[(i,12)] = round(math.sqrt(rms_s2n[(i,12)] / float(count_rms1)),6)
  rms_s2p[(i,12)] = round(math.sqrt(rms_s2p[(i,12)] / float(count_rms2)),6)
  # There's no S1n, S1p data for octupole table as of 02/26/2019.
  #rms_s1n[(i,12)] = round(math.sqrt(rms_s1n[(i,12)] / float(count_rms3)),6)
  #rms_s1p[(i,12)] = round(math.sqrt(rms_s1p[(i,12)] / float(count_rms4)),6)
  rms_be[(i,12)] = round(math.sqrt(rms_be[(i,12)] / float(count_rms5)),6)
  #print (mtN[i],"rms error of S2n:",rms_s2n[(i,12)],",nuclei count:",count_rms1)
  #print (mtN[i],"rms error of S2p:",rms_s2p[(i,12)],",nuclei count:",count_rms2)
  #print (mtN[i],"rms error of BE:",rms_be[(i,12)],",nuclei count:",count_rms5)
 # pikachu
 # Z,N = 2,8 wasn't calculated in the octupole table, Erik's range, causing some energies to be missing
# for Z in range(zs,zMax+1,2):
#   for N in range(ns,nMax+1,2):
#    if (N,Z,1,12) in ResS2p and (N,Z,23,12) not in ResS2p:
#        print (Z,N)

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
 return ResS2p, ResS2n, ResS1n, ResS1p, ResBE, newNucS1n, newNucS2n
 ###################################################################################################
 ###                                         4. PLOTTING                                         ###
 ###################################################################################################

 ###################################################################################################
 ### ResC: Dictionary of arrays for S2n/p_residual/raw plotting,                                 ###
 ### its key is structured as: (Z,i,j,k,l):                                                      ###
 ### Z: proton #;                                                                                ###
 ### i: theory masstable index;                                                                  ###
 ### j: experiment data index;                                                                   ###
 ### k: data type:: 0-S2n_res, 1-S2p_res; 2-S1n_res, 3-S1p_res; 4-BindingEnergy_res;             ###
 ### k: 5-S2n, 6-S2p, 7-BE                                                                       ###
 ### l= 0: neutron / mass / proton number sequence fox x-axis plotting                           ###
 ### l= 1: residual sequence for y-axis plotting                                                  ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n/p AS OF 5/31/18                                  ###
 ### Use ResC[(Z,****)] if you want to connect all isotopes as a chain in plotting later on      ###
 ### ResC[***].extend([N]) if you want to plot vs. N                                             ###
 ###################################################################################################
 #Dictionary of arrays for plotting purpose, each array key has 5 elements:
 #(proton #; theory masstable index; experiment dat index; data type 0: S2n, 1: S2p, 2: S1n; 0: neutron array, 1: residual array)


 # fFormat: figure format type when saving
 # pMod: plot mode "test" for testing, "global" or "local" for actual storing
 # pSave: save plot for 0-none; 1-S2n,BE_res,S2n_res; 2-S2p; 3-S2n & S2p; 4-new nuclei (AME2016-2003); 5-contour; 6-all
def plots(pSave = 0, fFormat = "pdf", pMod = "global"):

 S2n,S2p,S1n,S1p,BE,oct_nuclei = data_import()
 ResS2p, ResS2n, ResS1n, ResS1p, ResBE, newNucS1n, newNucS2n = residuals()

 ResC = {}
 ### All Isotonic Chain Goes Here ###
 ### All Isotonic Chain Goes Here ###
 ### All Isotonic Chain Goes Here ###
 zMax1 = 120
 # model 1~6, 14~16, 23~29 ; exp 12
 # oct_nuclei only works for the 7 Skyrme masstables
 for i in it.chain(range(1,7), [16], range(23,30)):
  for j in range(12,13):
   #FOR NOW WE ONLY PLOT EVEN-Z DATA
    for Z in range(2,zMax1+1,2):
     ResC[(Z,i,j,0,0)] = []; ResC[(Z,i,j,0,1)] = [] # S2n residual
     ResC[(Z,i,j,5,0)] = []; ResC[(Z,i,j,5,1)] = [] # S2n raw
     ResC[(Z,i,j,4,0)] = []; ResC[(Z,i,j,4,1)] = [] # BE residual
     ResC[(Z,i,j,7,0)] = []; ResC[(Z,i,j,7,1)] = [] # BE raw
     #FOR NOW WE ONLY PLOT EVEN-N DATA
     # Apply dripline cutoff
     n_min = dripline[(Z,i,0)] ; n_max = dripline[(Z,i,1)]
     # No dripline cutoff
     #n_min =2, n_max = nMax1
     for N in range(n_min,n_max+1,2):
      # S2n residual vs N
      if (N,Z,i) in oct_nuclei: flag = oct_nuclei[(N,Z,i)]
      else: flag = True
      if (N,Z,i,j) in ResS2n and flag:
       ResC[(Z,i,j,0,0)].extend([N])
       ResC[(Z,i,j,0,1)].extend([ResS2n[(N,Z,i,j)]])
      # S2n vs N, here j is a dummy index, we use 12-AME2016 just for convenience later on
      if (N,Z,i) in S2n and flag:
       ResC[(Z,i,j,5,0)].extend([N])
       ResC[(Z,i,j,5,1)].extend([S2n[(N,Z,i)]])
      # BE residual vs N
      if (N,Z,i,j) in ResBE and flag:
       ResC[(Z,i,j,4,0)].extend([N])
       ResC[(Z,i,j,4,1)].extend([ResBE[(N,Z,i,j)]])
      # BE vs N, here j is a dummy index, we use 12-AME2016 just for convenience later on
      if (N,Z,i) in BE and flag:
       ResC[(Z,i,j,7,0)].extend([N])
       ResC[(Z,i,j,7,1)].extend([BE[(N,Z,i)]])

 #(neutron #; theory masstable index; experiment dat index; data type 0: S2n, 1: S2p; 0: neutron sequence, 1: residual sequence)

 ### All Isotonic Chain Goes Here ###
 ### All Isotonic Chain Goes Here ###
 ### All Isotonic Chain Goes Here ###

 #Data storage for residual of S2p vs Proton # plot
 # model 1~6, 14~16, 23~29 ; exp 12
 for i in it.chain(range(1,7), range(14,17), range(23,30)):
  for j in range(12,13):
   #FOR NOW WE ONLY PLOT EVEN-N DATA
    for N in range(2,nMax1+1,2):
     ResC[(N,i,j,1,0)] = []; ResC[(N,i,j,1,1)] = [] # S2p residual
     ResC[(N,i,j,6,0)] = []; ResC[(N,i,j,6,1)] = [] # S2p raw
     #FOR NOW WE ONLY PLOT EVEN-Z DATA
     for Z in range(2,zMax1+1,2):
      if (N,Z,i,j) in ResS2p :
       ResC[(N,i,j,1,0)].extend([Z])
       ResC[(N,i,j,1,1)].extend([ResS2p[(N,Z,i,j)]])
      # S2p vs Z, here j is a dummy index, we use 12-AME2016 just for convenience later on
      if (N,Z,i) in S2p:
       ResC[(N,i,j,6,0)].extend([Z])
       ResC[(N,i,j,6,1)].extend([S2p[(N,Z,i)]])

 ### All Isobaric Chain Goes Here ###
 ### All Isobaric Chain Goes Here ###
 ### All Isobaric Chain Goes Here ###


 ###################################################################################################
 ###                                     4.0. PLOT OPTIONS                                       ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

 # 7-2d histogram
 scatter_BE = 1; compare_old = 0
 # zMax,nMax is largest proton/neutron # in provided experimental value
 # Range of Z to plot (for S2n)
 zS = 2; zE = zMax1           #do not exceed range of zS=2 to zE=zMax
 # Range of N to plot (for S2p)
 nS = 2; nE = nMax1          #do not exceed range of nS=2 to nE=nMax
 # Contour plotting, gaussian smoothing or none
 gs_smooth = False
 # gs_range is how far neighboring even-even neutrons and protons to smooth with
 gs_range = 4
 # gs_sigma defaults to 2, defines how spread the gaussian smoothing should be
 gs_sigma = 4; gs_lambda = 0.5 / (gs_sigma**2)
 # Contour color scaling separation energy cutoff:
 c_min = -1.5; c_max = 1.5
 ###################################################################################################
 ###                                     4.1. PLOTTING STARTS                                    ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 # Used to calculated chi-square between two sets of data
 s2n_rms = {}
 ###################################################################################################
 ### ResC                                                                                        ###
 ### k: data type:: 0-S2n_res, 1-S2p_res; 2-S1n_res, 3-S1p_res; 4-BindingEnergy_res;             ###
 ### k: 5-S2n, 6-S2p, 7-BE                                                                       ###
 ###################################################################################################
 # model 1~6, 16, 23~29 ; exp 12
 if pSave == 1 or pSave == 3 or pSave == 6 or pSave == 2:
     for i in it.chain(range(1,7), [16]):
      cor = {1:23,2:24,3:25,4:26,5:27,6:28,16:29}
      for j in range(12,13):
       for k in [5]:
        plt.clf()
        #String for figure name and directory when saved
        fN = "plots/"+pMod+"/" + mtN[j]+"/"
 ###################################################################################################
 ###                               4.a. S2N residual PLOTTING                                    ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 ### ResC                                                                                        ###
 ### k: data type:: 0-S2n_res, 1-S2p_res; 2-S1n_res, 3-S1p_res; 4-BindingEnergy_res;             ###
 ### k: 5-S2n, 6-S2p, 7-BE                                                                       ###
 ###################################################################################################
        #k = 0, plotting S2n data:
        if ( k==0 or k == 5 or k == 4 and (pSave == 1 or pSave == 3 or pSave == 6) ):
         s2n_rms[(i,j,k)] = 0; s2n_rms[(cor[i],j,k)] = 0
         plt.xlabel("Neutron Number",fontsize=14)
         if k == 0:
            fN += "S2n/Residual_s2n_" + mtN[i] + "_" + mtN[j]
         if k == 5:
            if not compare_old: fN += "S2n/S2n_oct_only_" + mtN[i]
            else: fN += "S2n/S2n" + mtN[i]
         if k == 4:
            fN += "BE/BE_" + mtN[i]
         if compare_old:
            fN += "_compare"
         fN += "." + fFormat
         count1 = 0; count2=0      #counts how many nuclei data point is in the plot
         leg_c = 0
         for Z in range(zS,zE+1,2):
          # Connect isotopic chain
          if k == 0 or k == 5:
          # New octupole masstable
            plt.plot(np.array(ResC[(Z,cor[i],j,k,0)]), np.array(ResC[Z,cor[i],j,k,1]),
                'b-'if Z%4 or compare_old else 'r-',lw=0.2,label='octupole' if leg_c == 0 else "",
                marker='.' if k!=5 else None,ms=1 if k!=5 else None)
          # Old masstable
            if compare_old:
                plt.plot(np.array(ResC[(Z,i,j,k,0)]), np.array(ResC[Z,i,j,k,1]),'r-',lw=0.2,
                label='old' if leg_c == 0 else "",marker='.' if k!=5 else None,ms=1 if k!=5 else None)
           #scatter plot for BE residual only
          if k == 4:
            if scatter_BE:
          # Scatter plot
          # Old masstable
              if compare_old:
                plt.scatter(np.array(ResC[(Z,i,j,k,0)]), np.array(ResC[Z,i,j,k,1]),c='r',s=1
                ,label='old' if leg_c == 0 else "" )
              # New octupole masstable
              plt.scatter(np.array(ResC[(Z,cor[i],j,k,0)]), np.array(ResC[Z,cor[i],j,k,1]),c='b',s=1
              ,label='octupole' if leg_c == 0 else "")
            elif not scatter_BE:
                    # New octupole masstable
                plt.plot(np.array(ResC[(Z,cor[i],j,k,0)]), np.array(ResC[Z,cor[i],j,k,1]),
                    'b-'if Z%4 or compare_old else 'r-',lw=0.2,label='octupole' if leg_c == 0 else "",
                    marker='.',ms=1)
              # Old masstable
                if compare_old:
                    plt.plot(np.array(ResC[(Z,i,j,k,0)]), np.array(ResC[Z,i,j,k,1]),'r-',lw=0.2,
                    label='old' if leg_c == 0 else "",marker='.',ms=1,)


          #,marker='.',ms=1,label='octupole' if leg_c == 0 else "")
          leg_c = 1
          count1 += len(np.array(ResC[Z,i,j,k,0]))
          count2 += len(np.array(ResC[Z,cor[i],j,k,0]))
          for a in range(0,min(len(np.array(ResC[Z,i,j,k,1])),len(np.array(ResC[Z,cor[i],j,k,1])))):
           s2n_rms[(i,j,k)] += (ResC[Z,i,j,k,1][a])**2
           s2n_rms[(cor[i],j,k)] += (ResC[Z,cor[i],j,k,1][a])**2
         #Normalizing chi-square with nuclei count
         s2n_rms[(i,j,k)] = math.sqrt(s2n_rms[(i,j,k)]/count1)
         s2n_rms[(cor[i],j,k)] = math.sqrt(s2n_rms[(cor[i],j,k)]/count2)
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
          # S2n residual
          if k == 0:
              xx = np.arange(0,161,20); yy = np.arange(-10,11,2)
              plt.title("$S_{2n}$ residual:  "+mtN[i]+" vs "+ mtN[j])
              plt.ylabel("$S_{2n,th}$ - $S_{2n,exp}$ (MeV)",fontsize=14)
              #On the top right corner, print data points and chi-square
              plt.text(115,8.5,"nuclei count: "+str(count1) )
              plt.text(115,7.5,"rms_old: " + str(round( s2n_rms[(i,j,k)],6 )) )
              plt.text(115,6.5,"rms_new: " + str(round( s2n_rms[(cor[i],j,k)],6 )) )
              plt.text(115,5.5,"improve: " + str(round( s2n_rms[(i,j,k)]-s2n_rms[(cor[i],j,k)],6 )) )
              plt.legend(loc='lower right')
          # S2n
          elif k == 5:
              xx = np.arange(0,301,20); yy = np.arange(0,61,10)
              axes = plt.gca(); axes.set_ylim([-5,60])
              plt.ylabel("$S_{2n}$ (MeV)",fontsize=14)
              plt.axvline(x=184,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.axvline(x=258,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.text(220,55,"nuclei count: "+str(count2) )
          # BE residual
          elif k == 4:
              xx = np.arange(0,161,20); yy = np.arange(-20,21,5)
              axes = plt.gca(); axes.set_ylim([-20,20])
              plt.text(110,16,"nuclei count: "+str(count1) )
              plt.text(110,12,"rms_new: " + str(round( s2n_rms[(cor[i],j,k)],6 )) )

              if compare_old:
                plt.text(110,14,"rms_old: " + str(round( s2n_rms[(i,j,k)],6 )) )
                plt.text(110,10,"improve: " + str(round( s2n_rms[(i,j,k)]-s2n_rms[(cor[i],j,k)],6 )) )
              plt.ylabel("Binding Energy residual(MeV)",fontsize=14)

          if k == 5 or k == 4:
              plt.axvline(x=20,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.axvline(x=28,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.axvline(x=50,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.axvline(x=82,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.axvline(x=126,linewidth=0.1,linestyle=':',color='k',alpha=0.7)
              plt.title(mtN[i])
              plt.legend(loc='upper left')
              #On the top right corner, print data points and chi-square


         plt.xticks(xx); plt.yticks(yy)
         #draw horizontal dashed line of 0 keV
         plt.axhline(y=0,linestyle='--',color='k')
         #save plot as fN
         plt.savefig(fN,format=fFormat)
         plt.close()
 ###################################################################################################
 ###                               4.a. S2P residual PLOTTING                                    ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1p AS OF 5/31/18                                    ###
 ###################################################################################################
        #k = 1, plotting S2p data:
        if ( k==1 or k == 6 and (pSave == 2 or pSave == 3 or pSave == 6) ):
         s2n_rms[(i,j,k)] = 0
         plt.xlabel("Proton Number",fontsize=14)
         if k == 1:
            fN = fN + "S2p/" + mtN[i] + "_" + mtN[j] + "_S2p_residual" + "." + fFormat
         if k == 6:
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
         #draw horizontal dashed line of 0 keV
         plt.axhline(y=0,linestyle='--',color='k')
         #save plot as fN
         plt.savefig(fN, format=fFormat)
         plt.close()
 ###################################################################################################
 ###                                 4.c. NEW NUCLEI POINTS PLOT                                 ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 if pSave == 4 or pSave == 6:
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
      if k == 0 :
       plt.plot( np.array(oldNS2n), np.array(oldZS2n), 'o', c='#33FBFF', ms = 3, label = "AME2003" )
       plt.plot( np.array(newNS2n), np.array(newZS2n), 'o', c='#087DF7', ms = 3, label = "AME2016")
       plt.plot( np.array(newerNS2n), np.array(newerZS2n), '*', c='#F70819', ms = 4, label = "2017" )
       fN = "plots/"+pMod+"/newS2n_nuclei.pdf"
       plt.title(r"Experimental S$_{2n}$ data")
       plt.legend(loc = 4)
       plt.savefig(fN, format=fFormat)
       plt.close()
      elif k == 1 :
       plt.plot( np.array(newNS2p), np.array(newZS2p), 'o', c='#C94631', ms = 5 )
       plt.title("AME2016 vs AME 2003 new S2p nuclei")
       fN = "plots/"+pMod+"/newS2p_nuclei."+fFormat
       plt.savefig(fN, format=fFormat)
       plt.close()
      #Plot both new S2n/S2p nuclei on one chart
      elif k == 2 :
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
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

 # ContourD: dictionary for contour plotting
 # key is: (theory masstable index 'i'; exp. masstable index 'j'; datatype k; array content: l)
 # For l=0: array of neutron numbers, later to be used on x-axis; l=1: array of proton numbers, later to be used on y-axis, l=3: array of data for bin coloring.
 # k : 0-S2n, 1-S2p, 2-beta2, 3-beta3
 contourD = {}
 gsN = 2; gsZ = 2;
 # model 1~10, 14,15,17; exp 11~13
# for i in it.chain(range(1,11), range(14,17)):
#  for j in range(11,14):
#   for k in range(0,4):
 for i in range(23,30):
  for j in range(12,13):
   for k in range(0,4):
    #for each i,j, initialize a new array for storing neutron/proton/data sequence for scipy contour plotting
    contourD[(i,j,k,0)] = []
    contourD[(i,j,k,1)] = []
    contourD[(i,j,k,2)] = []
    for Z in range(2,zMax1+1,2):
     for N in range(2,nMax1+1,2):
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
       # octupole table deformation
       if k == 2 or k == 3:
        if ( (N,Z,i) in deform ):
         contourD[(i,j,k,0)].extend([N])
         contourD[(i,j,k,1)].extend([Z])
         contourD[(i,j,k,2)].extend([deform[(N,Z,i)][k-2]]) #deform[(N,Z,i)] = (beta2,beta3)

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
 ###                          4.d.1 2D CONTOUR PLOT-6 Figures Combined                           ###
 ### THIS SECTION HAS NOT BEEN MODIFIED FOR S1n AS OF 5/31/18                                    ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 plt.clf()
 # deform[(N,Z,i)] = (float(data1), float(data2))
 # model 1~10, 14,15; exp 11~13
 # Models subplot UNEDF1, DDME2, FRDM2012, SLY4, DDPC1, HFB24, left to right, up to down
 if pSave == 5 or pSave == 6:
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
        if k==0 :
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
 ###                      4.d.2 2D HISTOGRAM PLOT-INDIVIDUAL MASSTABLE                           ###
 ### CREATED 2/27/2019 MAXWELL CAO                                                               ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 # plot beta2,beta3, j=12 dummy index do not change unless you want to change how contourD was created
 plt.clf()
 plt.close("all")
 if pSave == 7 or pSave == 6:
    for i in range(23,30):
        for k in [2,3]:
            x = np.array(contourD[(i,j,k,0)])
            y = np.array(contourD[(i,j,k,1)])
            fig, ax = plt.subplots(1,1,figsize=(6,2.3))
            ax2 = fig.add_axes([0.86, 0.1, 0.02, 0.7])
            if k == 2:
                cmap = plt.cm.seismic
                #cmaplist = [cmap(i) for i in range(cmap.N)]
                #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
                z = np.array(contourD[(i,j,k,2)])
                fN = "plots/2d_histogram/beta2_" + mtN[i]
                # define the bins
                bounds = np.linspace(-0.55,0.55,12)
                norm = BoundaryNorm(bounds, cmap.N)
                if refined_gs:
                    ax.set_title(mtN[i] + r" even-even $\beta_2$ refined")
                    fN += "_refined"
                else:
                    ax.set_title(mtN[i] + r" even-even $\beta_2$")
                ax2.set_title(r"   $\beta_2$", size=10)#, labelpad=0.5)
            elif k == 3:
                cmap = plt.cm.binary
                z = np.abs(np.array(contourD[(i,j,k,2)]))
                fN = "plots/2d_histogram/beta3_" + mtN[i]
                # define the bins
                cmap.set_under("w")
                bounds = np.linspace(0.0,0.5,21)
                norm = BoundaryNorm(bounds, cmap.N)
                if refined_gs:
                    ax.set_title(mtN[i] + r" even-even $\beta_3$ refined")
                    fN += "_refined"
                else:
                    ax.set_title(mtN[i] + r" even-even $\beta_3$")
                ax2.set_title(r"   $\beta_3$", size=10)#, labelpad=0.5)
            # colorbar
            cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds
            ,boundaries=bounds)
            # scatter plot
            ax.scatter(x, y, c=z, marker="s",linewidths=0.1,s=4,cmap=cmap,norm=norm,edgecolor='k')
            # create a second axes for the colorbar
            cb.ax.tick_params(labelsize=6)
            ax.set_yticks(range(0,121,20))
            ax.set_xticks(range(0,301,20))
            ax.tick_params(labelsize=8)
            plt.gcf().subplots_adjust(left=0.07,right=0.85)
            fN += "." +fFormat
            plt.savefig(fN,format=fFormat)




 ###################################################################################################
 ###                             5.a. NUMERIC DATA OUTPUT, residual                              ###
 ###################################################################################################
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

def write_output(saveNum=0, odevity=0, tFormat='csv'):
 S2n,S2p,S1n,S1p,BE,oct_nuclei,zMax,nMax = data_import()
 ResS2p, ResS2n, ResS1n, ResS1p, ResBE, newNucS1n, newNucS2n = residuals()

 if ( saveNum == 5 ):
  outputLabel = "Z".ljust(5) + "N".ljust(5) +  "SKM*(MeV)".ljust(CC) +  "SKP(MeV)".ljust(CC) +  "SLY4(MeV)".ljust(CC) +  "SVMIN(MeV)".ljust(CC) +  "UNEDF0(MeV)".ljust(CC) +  "UNEDF1(MeV)".ljust(CC) +  "UNEDF2(MeV)".ljust(CC) +  "DDME2(MeV)".ljust(CC) +  "DDMEd(MeV)".ljust(CC) +  "DDPC1(MeV)".ljust(CC) +  "NL3*(MeV)".ljust(CC) + "FDRM2012(MeV)".ljust(CC) +  "HFB24(MeV)".ljust(CC)
  outputStr0 = ""; outputStr1 = ""
  #count0, count1:scalars, count how many theories have residual for a specific nuclei
  #if all 10 theories don't have residual, count=0, then don't write that row of Z,N
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
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 ##
 ## NEUTRON SEPARATIONS
 ##
 seq = ("Z","N","AME2003","AME2003_sd","AME2016","AME2016_sd","JYFLTRAP2017","JYFLTRAP2017_sd","TRIUMF2018","TRIUMF2018_sd",                  "RIKEN2018","RIKEN2018_sd","new_other","new_other_sd","SKM*","SKP","SLY4","SVMIN","UNEDF0","UNEDF1","UNEDF2","FRDM2012","HFB24")
 outputLabel = ",".join(seq)
 #S1n
 if ( saveNum == 1):
  outputStr = ""
  #count0, count1:scalars, count how many theories have residual for a specific nuclei
  #if all theories don't have residual, count=0, then don't write that row of Z,N
  count = 0
  # even-even
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S1n_2018_even_Z_even_N."+tFormat, "w")
  # even Z, odd N
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S1n_2018_even_Z_odd_N."+tFormat, "w")
  # odd Z, even N
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S1n_2018_odd_Z_even_N."+tFormat, "w")
  # odd Z, odd N
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S1n_2018_odd_Z_odd_N."+tFormat, "w")
  # All S1n available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/data_S1n_2018_all."+tFormat, "w")
  output.write(outputLabel)

  for Z in range(zS,zMax1+1,zOdd):
   for N in range(nS,nMax1+1,nOdd):
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF2018, 18-RIKEN2018, 19-new_other
    for i in it.chain(range(11,14), [17,18,19]):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S1nErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
      #if i == 19: print (Z,N,"\t S1n:\t",S1n[(N,Z,i)])
     else:
      outputStr = outputStr + "*" + "," + "*" + ","
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
     else:
      if Z*1.0 > 0.68*N:
       outputStr = outputStr + "**" + "," #.ljust(CC)
      elif Z*1.0 <= 0.68*N:
       outputStr = outputStr + "*" + ","
    for i in range(15,16):
     if ( (N,Z,i) in S1n ):
      outputStr = outputStr + str(round(S1n[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
     else:
      if Z*1.0 > 0.68*N:
       outputStr = outputStr + "**" #.ljust(CC)
      elif Z*1.0 <= 0.68*N:
       outputStr = outputStr + "*"
    #only write to output if at least 1 set has S1n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr+"\n")
  output.close()

 #S2n
 if ( saveNum == 2):
  outputStr = ""
  #count0, count1:scalars, count how many theories have residual for a specific nuclei
  #if all 10 theories don't have residual, count=0, then don't write that row of Z,N
  count = 0
  #remember to change name of output file according to AME data year
  # even-even
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S2n_2018_even_Z_even_N."+tFormat, "w")
  # even Z, odd N
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S2n_2018_even_Z_odd_N."+tFormat, "w")
  # odd Z, even N
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S2n_2018_odd_Z_even_N."+tFormat, "w")
  # odd Z, odd N
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S2n_2018_odd_Z_odd_N."+tFormat, "w")
  # All S2n available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/data_S2n_2018_all."+tFormat, "w")
  output.write(outputLabel)

  points = {}
  for i in range(1,20):
    points[(mtN[i])] = 0

  for Z in range(zS,zMax1+1,zOdd):
   for N in range(nS,nMax1+1,nOdd):
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0;
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF2018, 18-RIKEN2018, 19-new_other
    for i in it.chain(range(11,14), [17,18,19]):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S2nErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
      #if i == 19: print (Z,N,"\t S2n:\t",S2n[(N,Z,i)])
     else:
      outputStr = outputStr + "*" + "," + "*" + ","
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      if Z*1.0 > 0.68*N:
       outputStr = outputStr + "**" + "," #.ljust(CC)
      elif Z*1.0 <= 0.68*N:
       outputStr = outputStr + "*" + ","
    for i in range(15,16):
     if ( (N,Z,i) in S2n ):
      outputStr = outputStr + str(round(S2n[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      if Z*1.0 > 0.68*N:
       outputStr = outputStr + "**"
      elif Z*1.0 <= 0.68*N:
       outputStr = outputStr + "*"
    #only write to output if at least 1 set has S2n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr+"\n")
  output.close()
 ##
 ## PROTON SEPARATIONS
 ##
 #S1p
 if ( saveNum == 3):
  outputStr = ""
  #count0, count1:scalars, count how many theories have residual for a specific nuclei
  #if all theories don't have residual, count=0, then don't write that row of Z,N
  count = 0
  # even-even
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S1p_2018_even_Z_even_N."+tFormat, "w")
  # even Z, odd N
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S1p_2018_even_Z_odd_N."+tFormat, "w")
  # odd Z, even N
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S1p_2018_odd_Z_even_N."+tFormat, "w")
  # odd Z, odd N
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S1p_2018_odd_Z_odd_N."+tFormat, "w")
  # All S1p available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/data_S1p_2018_all."+tFormat, "w")
  output.write(outputLabel)

  for Z in range(zS,zMax1+1,zOdd):
   for N in range(nS,nMax1+1,nOdd):
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF2018, 18-RIKEN2018, 19-new_other
    for i in it.chain(range(11,14), [17,18,19]):
     if ( (N,Z,i) in S1p ):
      outputStr = outputStr + str(round(S1p[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S1pErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
      #if i == 19: print (Z,N,"\t S1p:\t",S1p[(N,Z,i)])
     else:
      outputStr = outputStr + "*" + "," + "*" + ","
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S1p ):
      outputStr = outputStr + str(round(S1p[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
     else:
      if Z*1.0 < 0.68*N:
       outputStr = outputStr + "**" + "," #.ljust(CC)
      elif Z*1.0 >= 0.68*N:
       outputStr = outputStr + "*" + ","

    for i in range(15,16):
     if ( (N,Z,i) in S1p ):
      outputStr = outputStr + str(round(S1p[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
     else:
      if Z*1.0 < 0.68*N:
       outputStr = outputStr + "**"
      elif Z*1.0 >= 0.68*N:
       outputStr = outputStr + "*"
    #only write to output if at least 1 set has S1n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr+"\n")
  output.close()

 #S2p
 if ( saveNum == 4):
  outputStr = ""
  #count0, count1:scalars, count how many theories have residual for a specific nuclei
  #if all 10 theories don't have residual, count=0, then don't write that row of Z,N
  count = 0
  #remember to change name of output file according to AME data year
  # even-even
  if odevity == 0 :
      zS = 2; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S2p_2018_even_Z_even_N."+tFormat, "w")
  # even Z, odd N
  elif odevity == 1 :
      zS = 2; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S2p_2018_even_Z_odd_N."+tFormat, "w")
  # odd Z, even N
  elif odevity == 2 :
      zS = 3; nS = 2; zOdd = 2; nOdd = 2;
      output = open("data/data_S2p_2018_odd_Z_even_N."+tFormat, "w")
  # odd Z, odd N
  elif odevity == 3 :
      zS = 3; nS = 3; zOdd = 2; nOdd = 2;
      output = open("data/data_S2p_2018_odd_Z_odd_N."+tFormat, "w")
  # All S2p available
  elif odevity == 4 :
      zS = 2; nS = 2; zOdd = 1; nOdd = 1;
      output = open("data/data_S2p_2018_all."+tFormat, "w")
  output.write(outputLabel)

  points = {}
  zMax1 = 120; nMax1 = 300
  for i in range(1,20):
    points[(mtN[i])] = 0

  for Z in range(zS,zMax1+1,zOdd):
   for N in range(nS,nMax1+1,nOdd):
    outputStr = "\n" + str(Z) + "," + str(N) + ","
    count = 0;
    # experimental data, 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-TRIUMF2018, 18-RIKEN2018, 19-new_other
    for i in it.chain(range(11,14), [17,18,19]):
     if ( (N,Z,i) in S2p ):
      outputStr = outputStr + str(round(S2p[(N,Z,i)]+0.00000001,6) ) + "," + str(round(S2pErr[(N,Z,i)]+0.00000001,6) ) + ","
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
      #if i == 19: print (Z,N,"\t S2p:\t",S2p[(N,Z,i)])
     else:
      outputStr = outputStr + "*" + "," + "*" + ","

    # Theory model data
    for i in it.chain(range(1,7), [16,14]):
     if ( (N,Z,i) in S2p ):
      outputStr = outputStr + str(round(S2p[(N,Z,i)]+0.000000001,6) ) + "," #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      if Z*1.0 < 0.68*N:
       outputStr = outputStr + "**" + "," #.ljust(CC)
      elif Z*1.0 >= 0.68*N:
       outputStr = outputStr + "*" + ","

    for i in range(15,16):
     if ( (N,Z,i) in S2p ):
      outputStr = outputStr + str(round(S2p[(N,Z,i)]+0.000000001,6) ) #.ljust(CC)
      count = count + 1
      points[mtN[i]] = points[mtN[i]] + 1
     else:
      if Z*1.0 < 0.68*N:
       outputStr = outputStr + "**"
      elif Z*1.0 >= 0.68*N:
       outputStr = outputStr + "*"
    #only write to output if at least 1 set has S2n data for (N,Z)
    if (count != 0 ):
     output.write(outputStr+"\n")
  output.close()



 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
 #!!                                         Data Indexing                                       !!#
 #!! THEORY_OLD:                                                                                 !!#
 #!! 1-SKMS, 2-SKP, 3-SLY4, 4-SVMIN, 5-UNEDF0, 6-UNEDF1, 7-DDME2, 8-DDMED, 9-DDPC1, 10-NL3S,     !!#
 #!! 14-FRDM2012, 15-HFB24, 16-UNEDF2                                                            !!#
 #!!                                                                                             !!#
 #!! THEORY_OCTUPOLE_TABLE:                                                                      !!#
 #!! 23-SKMS, 24-SKP, 25-SLY4, 26-SVMIN, 27-UNEDF0, 28-UNEDF1, 29-UNEDF2                         !!#
 #!!                                                                                             !!#
 #!! EXPERIMENTS:                                                                                !!#
 #!! 11-AME2003, 12-AME2016, 13-JYFLTRAP2017, 17-AME2003_reserved, 18-AME2016_reserved           !!#
 #!! 19-TRIUMF2018, 20-RIKEN2018, 21-NEW_OTHER, 22-AME2016_MassExcess                            !!#
 #!!                                                                                             !!#
 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# 05/02/2019
# generate .csv file with columns: Z,N, S2p, S1p
# constraints: even Z, S2p<0; these are used to determine model mixing weights for two proton emitters

# old theory tables only
# S#*, BE, oct_nuclei keys: (N,Z,i)
# ResS#* keys: (N,Z,i,j)
def twoProton_weight_nuclei():
    S2n,S2p,S1n,S1p,BE,oct_nuclei,zMax,nMax = data_import()
    output = open('2p_emitter_weight_nuclei.csv','w')
    outStr = "Z,N,Exp.,S2p,S1p\n"
    for Z in range(2,zMax+1,2):
        for N in range(2,nMax+1):
            # Use AME2016 and later experimental data
            line = ''
            for i in [11,12,13,19,20,21]:
                if (N,Z,i) in S2p:
                    if S2p[(N,Z,i)] <= 0:
                        line = [str(Z),str(N),mtN[i],str(S2p[(N,Z,i)])]
                        if (N,Z,i) in S1p:
                            line.append(str(S1p[(N,Z,i)]))
            if line:
                line = ",".join(line)
                outStr += line + "\n"
    # Manually add 67Kr
    outStr += "36,31,Others,-1.69"
    output.write(outStr)
    output.close()

# for each theory or experiment, generate set of tuples (Z,N), that stores nuclei which have all 4 separation
# energies positive, meaning they are stable
def qa_filter_nuclei():
    pass


 ###################################################################################################
 ###                           6. DRIPLINES & S2N AVAILABLE REGIONS                              ###
 ###################################################################################################
 # Label nuclei with S2n or S2p data in 2016 with

 ###################################################################################################
 ###                                      FUNCTION EXECUTION                                     ###
 ###################################################################################################

 ###################################################################################################
 ###                               5. NUMERIC DATA OUTPUT, OPTIONS                               ###
 ###################################################################################################
 # savenum switch: 0-No save; 1-S1n, 2-S2n, 3-S1p, 4-S2p, 5-residual
 # Do saveNum = 1 / odevity = 1,3 ; S1n for odd neutron
 # Do saveNum = 2 / odevity = 0,2 ; S2n for even neutron
 # Do saveNum = 3 / odevity = 2,3 ; S1p for odd proton
 # Do saveNum = 4 / odevity = 0,1 ; S2p for even proton

# uncomment to output csv file of separation energies
#exec_list = [(1,1),(1,3),(2,0),(2,2),(3,2),(3,3),(4,0),(4,1)]
#for (arg1,arg2) in exec_list: write_output(arg1,arg2,"csv")
write_output( saveNum = 4, odevity = 0, tFormat = "csv")
twoProton_weight_nuclei()
