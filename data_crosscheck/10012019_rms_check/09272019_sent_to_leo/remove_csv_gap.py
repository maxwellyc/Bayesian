import os
file_list = []
for f in os.listdir():
    if f.endswith(".csv"):
        file_list.append(f)

for f in file_list:
    lines = open(f,'r').readlines()
    output = open(f.split(".")[0]+"_no_gaps.csv",'w')
    outstr = ""
    for ind, line in enumerate(lines):
        if not ind%2:
            outstr += line
    output.write(outstr)
