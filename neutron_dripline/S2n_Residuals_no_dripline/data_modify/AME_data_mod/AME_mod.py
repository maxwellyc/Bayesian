#For calculation of odd values please refer to: http://massexplorer.frib.msu.edu/content/masstables/Odd_Values_from_Even_Data.pdf

def readFile(DataFileIn,DataFileOut):
 # These are the locations of Neutron #, Proton #, Binding Energy, Neutron pairing energy, Proton pairing energy in the RMF mass table files from: "Global performance of covariant energy density functionals: Ground state observables of even-even nuclei and the estimate of theoretical uncertainties", Physical Review C 89, 054320 (2014)
 aa = 0
 az = 2
 as2n = 3
 as2nErr = 4
 as2p = 5
 as2pErr = 6
 asQa = 7
 asQaErr = 8
 
 CC = 13 #define padding for each colomn
 f1 = open(str(DataFileIn))
 output = open(str(DataFileOut), "w")
 lines = f1.readlines()
 #BindE={}                          #Creates an empty dictionary for binding energy
 #S1p = {}                          #Creates an empty dictionary for 1 proton separation energy
 S2p = {}                          #Creates an empty dictionary for 2 proton separation energy
 S2pErr = {}
 S2n = {}                          #Creates an empty dictionary for 2 neutron separation energy
 S2nErr = {}
 #These flags are here in case when, eg. S2n, data is marked by *,
 #then the position of S2p data in the line is moved forward by 1,
 #because then the error of S2n that meant to take one displacement is missing
 flagN = 0
 flagP = 0
 flagQa = 0
 A = 0
 Qa = {}                           #Creates an empty dictionary for Q_alpha value (Q value for alpha decay)
 QaErr = {}
 # function that detects if arg. is a number
 def isNum(s):
  try:
    float(s)
    return True
  except ValueError:
    return False
 nMax = 2
 zMax = 2
 for line in lines:
   ss = line.split()
    #print(ss)
   try:
     # This if statement is to catch cases where the first column is nuclear element symbol instead of mass number
     # This is caused by the experimental data having a "0" in the first column when a new mass number occurs, thus
     # after mass number 100, the first element becomes "0100" instead of one step earlier "0 99" which was separable
     # Fixed 04/10/18
     if (ss[int(aa)].isalpha() and ss[int(aa)] != "A"):
      print (ss[int(aa)])
      print (A)
      A = A+1
      aThis = A
      print (A)
      print ("**********************")
      az = 1
      as2n = 2
      as2nErr = 3
      as2p = 4
      as2pErr = 5
      asQa = 6
      asQaErr = 7
     else:
      A = int(float(ss[int(aa)])+0.0001)      #Read mass number of line
      aThis = A
      az = 2
      as2n = 3
      as2nErr = 4
      as2p = 5
      as2pErr = 6
      asQa = 7
      asQaErr = 8
     Z = int(float(ss[int(az)])+0.0001)      #Read proton number of line
     zThis = Z
     N = A - Z
     nThis = N
     if (nThis > nMax):             #Record maximum neutron/proton number exsiting in data set
      nMax = nThis
     if (zThis > zMax):
      zMax = zThis
     #S2n
     data = str(ss[int(as2n)])
     if (isNum(data)):
      #See purpose of flags at the beginning
      flagN = 0
      S2n[(N,Z)] = round(float(data)/1000.0,6)
      #S2n Error
      data = str(ss[int(as2nErr)])
      if (isNum(data)):
       S2nErr[(N,Z)] = round(float(data)/1000.0,6)
      else:
       S2nErr[(N,Z)] = "*"
     else:
      S2n[(N,Z)] = "*"
      S2nErr[(N,Z)] = "*"
      if ("#" not in str(ss[int(as2n)]) ):
       flagN = 1
     #S2p
     data = str(ss[int(as2p)-flagN])
     if (isNum(data)):
      flagP = 0
      S2p[(N,Z)] = round(float(data)/1000.0,6)
      data = str(ss[int(as2pErr)-flagN])
      #S2p Error
      if (isNum(data)):
       S2pErr[(N,Z)] = round(float(data)/1000.0,6)
      else:
       S2pErr[(N,Z)] = "*"
     else:
      S2p[(N,Z)] = "*"
      S2pErr[(N,Z)] = "*"
      if ("#" not in str(ss[int(as2p)]) ):
       flagP = 1
     #Q_alpha
     data = str(ss[int(asQa)-flagN-flagP])
     if (isNum(data)):
      flagQa = 0
      Qa[(N,Z)] = round(float(data)/1000.0,6)
      data = str(ss[int(asQaErr)-flagN-flagP])
      #Q_alpha error
      if (isNum(data)):
       QaErr[(N,Z)] = round(float(data)/1000.0,6)
      else:
       QaErr[(N,Z)] = "*"
     else:
      Qa[(N,Z)] = "*"
      QaErr[(N,Z)] = "*"
      if ("#" not in str(ss[int(asQa)]) ):
       flagQa = 1
   except (ValueError, IndexError):        #N,Z, or, S2n/S2p are not numbers
       continue
 outputStr = "Z".ljust(5) + "N".ljust(5) +  "S1n(MeV)".ljust(CC+5) + "S1n_Err(MeV)".ljust(CC+5) +\
  "S1p(MeV)".ljust(CC+5) + "S1p_Error(MeV)".ljust(CC+5) + "Q4b(MeV)".ljust(CC+5) +\
  "Q4b_Error(MeV)".ljust(CC+5)
#   outputStr = "Z".ljust(5) + "N".ljust(5) +  "S2n(KeV)".ljust(CC+5) + "S2n_Err(MeV)".ljust(CC+5) +\
#  "S2p(MeV)".ljust(CC+5) + "S2p_Error(MeV)".ljust(CC+5) + "Qa(MeV)".ljust(CC+5) +\
#  "Qa_Error(MeV)".ljust(CC+5)
 output.write(outputStr)
 
 for Z in range(1,zMax+1):
  for N in range(1,nMax+1):
   # Check if there's any of these data available for this N,Z
   if ( (N,Z) in S2n or (N,Z) in S2p or (N,Z) in Qa ):
    outputStr = "\n" + str(Z).ljust(5) + str(N).ljust(5)
    #Write S2n
    if ( (N,Z) in S2n ):
     outputStr = outputStr + str(S2n[(N,Z)]).ljust(CC+5)
     if ( (N,Z) in S2nErr ):
      outputStr = outputStr + str(S2nErr[(N,Z)]).ljust(CC+5)
     else:
      outputStr = outputStr + "*".ljust(CC+5)
    else:
     outputStr = outputStr + "*".ljust(CC+5) + "*".ljust(CC+5)
    #Write S2p
    if ( (N,Z) in S2p ):
     outputStr = outputStr + str(S2p[(N,Z)]).ljust(CC+5)
     if ( (N,Z) in S2pErr ):
      outputStr = outputStr + str(S2pErr[(N,Z)]).ljust(CC+5)
     else:
      outputStr = outputStr + "*".ljust(CC+5)
    else:
     outputStr = outputStr + "*".ljust(CC+5) + "*".ljust(CC+5)
    #Write Q_alpha
    if ( (N,Z) in Qa ):
     outputStr = outputStr + str(Qa[(N,Z)]).ljust(CC+5)
     if ( (N,Z) in QaErr ):
      outputStr = outputStr + str(QaErr[(N,Z)]).ljust(CC+5)
     else:
      outputStr = outputStr + "*".ljust(CC+5)
    else:
     outputStr = outputStr + "*".ljust(CC+5) + "*".ljust(CC+5)
    if ( S2n[(N,Z)] != "*" or S2p[(N,Z)] != "*" or Qa[(N,Z)] != "*"):
     output.write(outputStr)
    


 #print(S2p)
 f1.close()
 output.close()
DataFileIn='AME_NoHeader/rct2.mas03NH.txt'
DataFileOut='2003AME_S1n.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

DataFileIn='AME_NoHeader/rct2-16NH.txt'
DataFileOut='2016AME_S1n.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

