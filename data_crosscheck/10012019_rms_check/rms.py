import os
import math

def isNum(x):
    try:
        float(x)
        return True
    except:
        return False

dir = "09272019_sent_to_leo"
theory_start_column = 14

edf_names,file_list,types = [], [], set()
exp_data, theo_data, rms = {}, {}, {}

# Leo's list of nuclei
leo_set = set()
lines = open("nuclei_rms.csv",'r').readlines()[1:]
for line in lines:
    ss = line.split(",")
    Z, N = int(ss[2]), int(ss[3])
    leo_set.add((Z,N))

for f in os.listdir(dir):
    if f.endswith(".csv"): file_list.append(f)

for f in file_list:

    type_start, type_end = f.find("S"), f.find(".csv")
    file_type = f[type_start:type_start+3]
    types.add(file_type)
    lines = open("./"+dir+"/"+f,'r').readlines()
    edf_names = lines.pop(0).split(",")[theory_start_column:]
    if "\n" in edf_names[-1]:
        edf_names[-1] = edf_names[-1][:-1]

    if file_type not in theo_data:
        theo_data[file_type] = {}
        # each 2nd level dictionary belongs to an edf
        for edf in edf_names:
            if edf == "UNEDF1-SO": continue
            theo_data[file_type][edf] = {}

    if file_type not in exp_data:
        exp_data[file_type] = {}

    # go through entire file and extract data
    for line in lines:
        ss = line.split(",")
        Z, N = int(ss[0]), int(ss[1])

        # use newest available data
        for i in range(2,14,2):
            val = ss[i]
            if isNum(val):
                exp_data[file_type][(Z,N)] = float(val)

        # Only record theory data if exp. exsists
        # if (Z,N) in exp_data[file_type]:
        # Only record theory data if in Leo's set of nuclei for RMS computation
        if (Z,N) in leo_set:
            for ind,edf in enumerate(edf_names):
                if edf == "UNEDF1-SO": continue
                val = ss[theory_start_column+ind]
                if isNum(val):
                    theo_data[file_type][edf][(Z,N)] = float(val)

#print (theo_data["S2p"])
# RMS error print out
for file_type in ["S2p"]:
    #print (file_type, "Experimental d.p. count:",len(exp_data[file_type]))
    for edf in edf_names:
        sums, count = 0, 0
        if edf == "UNEDF1-SO": continue
        for Z,N in theo_data[file_type][edf]:
            sums += (theo_data[file_type][edf][(Z,N)] - exp_data[file_type][(Z,N)])**2
            count += 1
        rms[edf] = round(math.sqrt(sums / count),2)
        print (edf+"    ", "\t", rms[edf], "\tcounts:",count)
    print ("\n")

# debug section:
# sums_exp, count_exp = 0,0
# for Z,N in sorted(exp_data["S2p"]):
#     print (Z,N,exp_data["S2p"][(Z,N)])
#     sums_exp += exp_data["S2p"][(Z,N)]
#     count_exp += 1
#
# print (sums_exp, count_exp)



# for Z,N in sorted(theo_data["S2p"]["SKMS"]):
#     print (Z,N,theo_data["S2p"]["SKMS"][(Z,N)])
