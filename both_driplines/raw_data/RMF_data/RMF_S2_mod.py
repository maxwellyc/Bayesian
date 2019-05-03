#For calculation of odd values please refer to: http://massexplorer.frib.msu.edu/content/masstables/Odd_Values_from_Even_Data.pdf

def readFile(DataFileIn,DataFileOut):
 # These are the locations of Neutron #, Proton #, Binding Energy, Neutron pairing energy, Proton pairing energy in the Relativistic Mean Field mass table files from: "Global performance of covariant energy density functionals: Ground state observables of even-even nuclei and the estimate of theoretical uncertainties", Physical Review C 89, 054320 (2014)
 an = 1
 az = 0
 ae = 3
 apn = 6
 apz = 7
 CC = 13 #define padding for each colomn
 f1 = open(str(DataFileIn))
 output = open(str(DataFileOut), "w")
 lines = f1.readlines()
 BindE={}                          #Creates an empty dictionary for binding energy
 PairN = {}                        #Creates an empty dictionary for neutron pairing gap
 PairZ = {}                        #Creates an empty dictionary for proton pairing gap
 S1p = {}                          #Creates an empty dictionary for 1 proton separation energy
 S2p = {}                          #Creates an empty dictionary for 2 proton separation energy
 S1n = {}                          #Creates an empty dictionary for 1 neutron separation energy
 S2n = {}                          #Creates an empty dictionary for 2 neutron separation energy
 Qa = {}                           #Creates an empty dictionary for Q_alpha value (Q value for alpha decay)
 # function that detects if arg. is a number
 outputStr = "Z".ljust(5) + "N".ljust(5) +  "Binding_E_(MeV)".ljust(CC+5) +  "S_1p_(MeV)".ljust(CC) +  "S_2p_(MeV)".ljust(CC) +  "S_1n_(MeV)".ljust(CC) +  "S_2n_(MeV)".ljust(CC) + "Q_alpha_(MeV)".ljust(CC)
 output.write(outputStr)
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
     N = int(float(ss[int(an)])+0.0001)      #Number of Neutrons
     nThis = N
     Z = int(float(ss[int(az)])+0.0001)      #Number of Protons
     zThis = Z
     if (nThis > nMax):
      nMax = nThis
     if (zThis > zMax):
      zMax = zThis
     BindE[(N,Z)] = float(ss[int(ae)])                      #Binding Energy Dict
     PairN[(N,Z)] = -float(ss[int(apn)])
     PairZ[(N,Z)] = -float(ss[int(apz)])
     # For each even-even data on file, the following computes binding energies of its neighbours as follows:
     # (Z-1,N-1) *     * (Z-1, N)
     #             \   |
     #              \  |
     #  (Z, N-1) * --- * (even Z,even N)
     #
     #Compute Odd-Z Even-N Binding energies
     if (BindE.has_key((N,Z)) and BindE.has_key((N,Z-2)) and PairZ.has_key((N,Z)) and PairZ.has_key((N,Z-2)) ):
        BindE[(N,Z-1)] = 0.5 * ( BindE[(N,Z)] + BindE[(N,Z-2)] + PairZ[(N,Z)] + PairZ[(N,Z-2)] )
     #Compute Even-Z Odd-N Binding energies
     if (BindE.has_key((N,Z)) and BindE.has_key((N-2,Z)) and PairN.has_key((N,Z)) and PairN.has_key((N-2,Z)) ):
        BindE[(N-1,Z)] = 0.5 * ( BindE[(N,Z)] + BindE[(N-2,Z)] + PairN[(N,Z)] + PairN[(N-2,Z)] )
     #Compute Odd-N Proton Pairing gaps, this is required for Odd-Z Odd-N Binding energies computation
     if (PairZ.has_key((N,Z)) and PairZ.has_key((N-2,Z))):
        PairZ[(N-1,Z)] = 0.5 * ( PairZ[(N,Z)] + PairZ[(N-2,Z)] )
     #Compute Odd-Z Odd-N Binding energies, this code works because the mass table is ordered in Proton number first, thus the complete Z-2 data is guaranteed to exist when the following lines are executed
     if ( BindE.has_key((N-1,Z)) and BindE.has_key((N-1,Z-2)) and PairZ.has_key((N-1,Z)) and PairZ.has_key((N-1,Z-2)) ):
        BindE[(N-1,Z-1)] = 0.5 * ( BindE[(N-1,Z)] + BindE[(N-1,Z-2)] + PairZ[(N-1,Z)] + PairZ[N-1,Z-2]  )
     #!!!
     #!!! Upon this point, all odd-odd, odd-even, even-odd Binding energies are computed and stored !!!
     #!!!
    except (ValueError, IndexError):        #N,Z, or, BE are not numbers
       continue
 for Z in range(2,zMax+1):
  for N in range(2,nMax+1):
   if ( BindE.has_key((N,Z)) ):
    outputStr = "\n" + str(Z).ljust(5) + str(N).ljust(5) +  str(BindE[(N,Z)]).ljust(CC+5)
    if ( BindE.has_key((N,Z-1)) ):
     S1p[(N,Z)] = BindE[(N,Z-1)] - BindE[(N,Z)]
     outputStr = outputStr + str( round( S1p[(N,Z)]+0.00000001,6 ) ).ljust(CC)
    else:
     outputStr = outputStr + "*".ljust(CC)
    if ( BindE.has_key((N,Z-2)) ):
     S2p[(N,Z)] = BindE[(N,Z-2)] - BindE[(N,Z)]
     outputStr = outputStr + str( round( S2p[(N,Z)]+0.00000001,6 ) ).ljust(CC)
    else:
     outputStr = outputStr + "*".ljust(CC)
    if ( BindE.has_key((N-1,Z)) ):
     S1n[(N,Z)] = BindE[(N-1,Z)] - BindE[(N,Z)]
     outputStr = outputStr + str( round( S1n[(N,Z)]+0.00000001,6 ) ).ljust(CC)
    else:
     outputStr = outputStr + "*".ljust(CC)
    if ( BindE.has_key((N-2,Z)) ):
     S2n[(N,Z)] = BindE[(N-2,Z)] - BindE[(N,Z)]
     outputStr = outputStr + str( round( S2n[(N,Z)]+0.00000001,6 ) ).ljust(CC)
    else:
     outputStr = outputStr + "*".ljust(CC)
    if ( BindE.has_key((N-2,Z-2)) ):
     Qa[(N,Z)] = 28.3 + BindE[(N,Z)] - BindE[(N-2,Z-2)]
     outputStr = outputStr + str( round( Qa[(N,Z)]+0.00000001,6 ) ).ljust(CC)
    else:
     outputStr = outputStr + "*".ljust(CC)
    output.write(outputStr)

 f1.close()
 output.close()


DataFileIn='RMFnoHeader/ddme2-tableNH.dat'
DataFileOut='RMFCompleteTable/ddme2-sep.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

DataFileIn='RMFnoHeader/ddmed-tableNH.dat'
DataFileOut='RMFCompleteTable/ddmed-sep.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

DataFileIn='RMFnoHeader/ddpc1-tableNH.dat'
DataFileOut='RMFCompleteTable/ddpc1-sep.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

DataFileIn='RMFnoHeader/nl3s-tableNH.dat'
DataFileOut='RMFCompleteTable/nl3s-sep.dat'
readFile(DataFileIn,DataFileOut)    #arg: (DataFileIn,DataFileOut)

