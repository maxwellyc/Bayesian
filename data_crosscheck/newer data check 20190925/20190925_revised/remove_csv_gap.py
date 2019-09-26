lines = open("data_yc_20190925_S2p_even_Z_odd_N"+".csv",'r').readlines()
output = open("data_eoyc_2019_S2p_even_Z_odd_N"+".csv",'w')
outstr = ""
for ind, line in enumerate(lines):
    if not ind%2:
        outstr += line
output.write(outstr)
