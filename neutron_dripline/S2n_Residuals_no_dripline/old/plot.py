import numpy as np
import matplotlib.pyplot as plt


plt.ylabel("$S_{2n,th}$ - $S_{2n,exp}$ (KeV)",fontsize=14)
plt.xlabel(r'Neutron numbers',fontsize=14)
#plt.suptitle('Ca isotopes S2n Theory-Exp.') #number of non-zero elements
#Arrays compared with 2003/2016 AME data starts with a*/b*
#Array indices in order of SKMS,SKP,SLY4,SV-MIN,UNEDF0,UNEDF1,DDME2,DDMED,DDPC1,NL3*

#Neutron numbers
xx = np.array([16,18,20,22,24,26,28,30,32,34,36])
ax = np.array([18,20,22,24,26,28,30,32])
bx = np.array([18,20,22,24,26,28,30,32,34])

#2003 AME S2N ERROR
aErr = np.array([40.26,4.55,0.14,0.29,2.28,4.11,8.32,698.68])

#2016 AME S2N ERROR
bErr = np.array([40,0.20,0.15,0.29,2.3,2.2,1.6,1.7,50])

#SKMS-2003 S2n (KeV)
a0 = np.array([-2695.719,-2796.323,2002.164,1128.890,1494.716,229.009,2217.519,2265.545])
#SKMS-2016 S2n (KeV)
b0 = np.array([-2683.229,-2797.153,2002.104,1128.870,1490.646,222.799,2209.659,526.785,784.033])

#SKP-2003 S2n (KeV)
a1 = np.array([-1005.290,-2150.178,1630.082,29.417,-272.843,-1510.999,2223.511,2102.263])
#SKP-2016 S2n (KeV)
b1 = np.array([-992.800,-2151.008,1630.022,29.397,-276.913,-1517.209,2215.651,363.503,1097.645])

#SLY4-2003 S2n (KeV)
a2 = np.array([-771.384,-1407.575,1271.315,-137.241,30.416,-1506.945,502.883,695.450])
#SLY4-2016 S2n (KeV)
b2 = np.array([-758.894,-1408.405,1271.255,-137.261,26.346,-1513.155,495.023,-1043.310,231.772])

#SVMIN-2003 S2n (KeV)
a3 = np.array([-1063.259,-1911.176,1541.675,-15.490,-235.458,-1635.828,1583.850,1285.304])
#SVMIN-2016 S2n (KeV)
b3 = np.array([-1050.769,-1912.006,1541.615,-15.510,-239.528,-1642.038,1575.990,-453.456,629.938])

#UNEDF0-2003 S2n (KeV)
a4 = np.array([-2241.925,-2680.420,1667.719,177.734,-219.947,-1621.811,1851.930,1524.726])
#UNEDF1-2016 S2n (KeV)
b4 = np.array([-2229.435,-2681.250,1667.659,177.714,-224.017,-1628.021,1844.070,-214.034,517.515])

#UNEDF1-2003 S2n (KeV)
a5 = np.array([-60.092,-931.366,1288.416,-177.190,-448.362,-1960.579,1154.349,1000.933])
#UNEDF1-2016 S2n (KeV)
b5 = np.array([-47.602,-932.196,1288.356,-177.210,-452.432,-1966.789,1146.489,-737.827,32.052])

#DDME2-2003 S2n (KeV)
a6 = np.array([587.510,865.310,721.570,-254.040,-550.230,-1332.690,-679.340,468.060])
#DDME2-2016 S2n (KeV)
b6 = np.array([600.000,864.480,721.510,-254.060,-554.300,-1338.900,-687.200,-1270.700,1182.000])

#DDMED-2003 S2n (KeV)
a7 = np.array([1270.510,1408.310,758.570,-757.040,-831.230,-1652.690,-1174.340,265.060])
#DDMED-2016 S2n (KeV)
b7 = np.array([1283.000,1407.480,758.510,-757.060,-835.300,-1658.900,-1182.200,-1473.700,1481.000])

#DDPC1-2003 S2n (KeV)
a8 = np.array([816.510,303.310,1407.570,-335.040,-262.230,-1289.690,-179.340,644.060])
#DDPC1-2016 S2n (KeV)
b8 = np.array([829.000,302.480,1407.510,-335.060,-266.300,-1295.900,-187.200,-1094.700,1487.000])

#NL3*-2003 S2n (KeV)
a9 = np.array([272.510,145.310,894.570,-208.040,-575.230,-1961.690,409.660,1297.060])
#NL3*-2016 S2n (KeV)
b9 = np.array([285.000,144.480,894.510,-208.060,-579.300,-1967.900,401.800,-441.700,2135.000])



#ax1 = plt.subplot()


#ax2 = ax1.twinx()
#Array indices in order of 0-SKMS,1-SKP,2-SLY4,3-SV-MIN,4-UNEDF0,5-UNEDF1,6-DDME2,7-DDMED,8-DDPC1,9-NL3*
plt.ylim([-3000,3000])
plt.yticks(np.arange(-3000,3001,500))
plt.title('NL3*')
plt.plot(ax,a9,'yo-',label='2003',ms=6)
plt.plot(bx,b9,'b^--',label='2016',ms=5)
plt.legend(numpoints=1,loc=4,fontsize=10)
plt.xticks(xx)
plt.axhline(y=0,linestyle='--',color='k')
plt.text(17,2500,"Ca", fontsize=20 )
#ax1.plot(x0,y1, 'b-')
#ax1.plot(x,y2, 'r-')
#legend = ax1.legend(loc='lower right', shadow=True,numpoints = 1)
#legend = plt.legend(loc='lower right', shadow=True,numpoints = 1)

plt.show()
