# read in AME data,
# check length of row array ss, must be 8, except 2 headers
# function to convert mass excess to binding energies
# double check with BE/A

# BE = Z * M_H + N * M_N - ( Mass excess + A * u)
def ME_BE(me,Z,N):
    # en.wikipedia.org/wiki/Nuclear_binding_energy
    # Example_values_deduced_from_experimentally_measured_atom_nuclide_masses
    # Below masses taken from above example value
    uAtom = 931494.028
    M_H = 938783.0802
    M_N = 939565.4133
    return round(Z*M_H + N*M_N - me - (Z+N) * uAtom,6)

def ame_mass_be(year):
    file = open('AME_NoHeader/AME'+str(year)+'_massNH.txt','r')
    lines = file.readlines()
    output = open(str(year)+'AME_mass.csv','w')
    outStr = "Z,N,BE(MeV),BE_sd(MeV),BE_BEA(MeV),diff(MeV)\n"
    count = 0
    BE_dict; BE_ERR_dict = {},{}
    for line in lines:
        # 2 rows of headers
        if count<2: count+=1 ; continue
        ss = line.split()
        if len(ss) != 8: print (f'Error: {ss}')
        A,Z,N = int(ss[2]),int(ss[1]),int(ss[0])
        me,meERR,be_a,be_aERR = ss[-4],ss[-3],ss[-2],ss[-1]
        if '#' in me: continue
        else:
            me,meERR,be_a,be_aERR = float(ss[-4]),float(ss[-3]),float(ss[-2]),float(ss[-1])
            # binding from BE/A, for verification only
            BE_BEA = be_a*A/1000.0
            # binding from mass excess, for actual BE value
            BE_ME = ME_BE(me,Z,N)/1000.0
            meERR /= 1000.0
            # difference between 2 methods
            diff = BE_ME - BE_BEA
            BE_dict[(Z,N)] = BE_ME
            BE_ERR_dict[(Z,N)] = meERR
        outStr += str(Z) + ',' + str(N) + ',' + str(BE_ME) + ',' + str(meERR) + ',' +\
         str(BE_BEA) + ',' + str(diff) + '\n'
    output.write(outStr)
    file.close(); output.close()

ame_mass_be(2003)
ame_mass_be(2016)



