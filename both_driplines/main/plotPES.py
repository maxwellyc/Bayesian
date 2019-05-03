#rewrite 4 potential energy surface .dat files (2x denser on each axis), into a list of lists, with [int Q2,int Q3,double HFBenergy] as format.
#16 cores * 2 Q2 value calc. per core * 13 Q3 value per Q2 value = 416 lines per .dat file, 416 * 4 = 1664 in total
#After getting list of lists, start plotting PES into contour.

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate



def PES(InFile0,InFile1,InFile2,InFile3):
 f0 = open(str(InFile0))
 f1 = open(str(InFile1))
 f2 = open(str(InFile2))
 f3 = open(str(InFile3))
 lines0 = f0.readlines()
 lines1 = f1.readlines()
 lines2 = f2.readlines()
 lines3 = f3.readlines()
 Q3list = [0]*1664
 Q2list = [0]*1664
 Elist = [0]*1664
 PESmat = []      #Creates an empty 3d matrix, set key to Q2,Q3

 for line in lines0:
     PESraw = line.split()
     Q2 = int(round(float(PESraw[5])+0.01,0))
     Q3 = int(round(float(PESraw[6])+0.01,0))
     HFB_E = float(PESraw[4])
     PESmat.append([Q2,Q3,HFB_E])
 
 for line in lines1:
     PESraw = line.split()
     Q2 = int(round(float(PESraw[5])+0.01,0))
     Q3 = int(round(float(PESraw[6])+0.01,0))
     HFB_E = float(PESraw[4])
     PESmat.append([Q2,Q3,HFB_E])

 for line in lines2:
     PESraw = line.split()
     Q2 = int(round(float(PESraw[5])+0.01,0))
     Q3 = int(round(float(PESraw[6])+0.01,0))
     HFB_E = float(PESraw[4])
     PESmat.append([Q2,Q3,HFB_E])

 for line in lines3:
     PESraw = line.split()
     Q2 = int(round(float(PESraw[5])+0.01,0))
     Q3 = int(round(float(PESraw[6])+0.01,0))
     HFB_E = float(PESraw[4])
     PESmat.append([Q2,Q3,HFB_E])

 # FIND Lowest energy
 PESmat_len = len(PESmat)
 HFB_min = PESmat[0]
 for i in range (1,PESmat_len):
     HFB_E = PESmat[i][2]
     if (HFB_E < HFB_min):
         HFB_min = HFB_E
 
 count = 0
 #print PESmat[150]
 # Adjusting HFB energy to a relative value compared to lowest energy. For plotting purposes
 for i in range (0,PESmat_len):
     HFB_E = PESmat[i][2] - HFB_min
     PESmat[i][2] = HFB_E
     if (HFB_E < 1.0):
       Q2list[count] = int(PESmat[i][0])
       Q3list[count] = int(PESmat[i][1])
       Elist[count] = float(PESmat[i][2])
       count = count + 1
 
 count0 = 0
 
 for i in range (0, PESmat_len):
     if (Elist[i] != 0):
       count0 = count0 + 1

 Q2list0 = [0]*count0
 Q3list0 = [0]*count0
 Elist0 = [0]*count0

 for i in range (0,count0):
     Q2list0[i] = Q2list[i]
     Q3list0[i] = Q3list[i]
     Elist0[i] = Elist[i]

 print Elist
 #print PESmat
#PRINT SMALLEST PES data point
 for i in range (0,PESmat_len):
     if (PESmat[i][2] == 0):
         print PESmat[i]


 # Generate data:
 x = np.array(Q2list0)
 y = np.array(Q3list0)
 z = np.array(Elist0)
# Set up a regular grid of interpolation points
 xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
 xi, yi = np.meshgrid(xi, yi)

# Interpolate
 rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
 zi = rbf(xi, yi)

 im = plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()],aspect='auto')
 plt.scatter(x, y, c=z, marker=".",color = 'black')
 cb = plt.colorbar()
 cb.set_label(r'$E_0$ - $E_{MIN}$ (MeV)')
 plt.annotate('Q2',xy=(5,5), xytext=(-50,-680), fontsize=15)
 plt.annotate('Q3',xy=(5,5), xytext=(-2000,2930), fontsize=15)
 plt.title('UNEDF2 Pu240')
 plt.show()


 #print PESmat


InFile0='Pu240_PES0.dat'
InFile1='Pu240_PES1.dat'
InFile2='Pu240_PES2.dat'
InFile3='Pu240_PES3.dat'

PES(InFile0,InFile1,InFile2,InFile3)
