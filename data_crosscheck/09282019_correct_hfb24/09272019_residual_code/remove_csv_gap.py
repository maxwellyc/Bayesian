import os
file_list = []
for f in os.listdir():
    if f.endswith(".csv"):
        file_list.append(f)

for f in file_list:
    lines = open(f,'r').readlines()
    outstr = ""
    for ind, line in enumerate(lines):
        if ind%2 or ind == 0:
            outstr += line
    output = open(f,'w')
    output.write(outstr)
    output.close()
