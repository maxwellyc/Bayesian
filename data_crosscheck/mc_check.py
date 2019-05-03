# read in AME data,
# check length of row array ss, must be 8, except 2 headers
# function to convert mass excess to binding energies
# double check with BE/A
import math
# BE = Z * M_H + N * M_N - ( Mass excess + A * u)
def ME_BE(me,Z,N):
    # +- 0.0057 keV
    uAtom = 931494.0954
    # +- 0.0058 keV
    M_P = 938272.0813
    # +- 0.0000031 keV
    M_e = 510.9989461
    M_H = M_P + M_e
    # +- 0.0058 keV
    M_N = 939565.4133

    return round(Z*M_H + N*M_N - me - (Z+N) * uAtom,12)

# returns dictionary with key (Z,N), value (BE, BE_ERR), energies are in MeV
def AME_me2be(year):
    file = open('raw_exp/AME'+str(year)+'_massNH.txt','r')
    lines = file.readlines()
    count = 0
    temp_d = {}
    # AME masses
    for line in lines:
        # 2 rows of headers
        if count<2: count+=1 ; continue
        ss = line.split()
        if len(ss) != 8: print (f'Error: {ss}')
        A,Z,N = int(ss[2]),int(ss[1]),int(ss[0])
        me,meERR = ss[-4],ss[-3]
        if '#' in me: continue
        else:
            me,meERR = float(ss[-4]),float(ss[-3])
            # binding from mass excess, for actual BE value
            BE_ME = round(ME_BE(me,Z,N)/1000.0,6)
            meERR = round(meERR/1000.0,6)
            temp_d[(Z,N)] = (BE_ME, meERR)
    file.close()
    return temp_d

def AME_qa(year):
    file = open('raw_exp/AME'+str(year)+'_qaNH.txt','r')
    lines = file.readlines()
    count = 0
    temp_d = {}
    # AME masses
    for line in lines:
        # 2 rows of headers
        if count<2: count+=1 ; continue
        ss = line.split()
        A,Z = int(ss[0]),int(ss[2])
        N = A - Z
        shift = 0
        # If S2n value is valid
        if ss[3] != '*':
            shift += 1
        # If S2p value is valid
        if ss[4+shift] != '*':
            shift += 1
        qa,qaERR = ss[5+shift],ss[6+shift]
        if qa != '*' and '#' not in qa:
            qa,qaERR = round(float(qa)/1000.0,6),round(float(qaERR)/1000.0,6)
            temp_d[(Z,N)] = (qa, qaERR)
        else: continue
    file.close()
    return temp_d

# returns dictionary with key (Z,N), value (BE, BE_ERR), energies are in MeV
def JYFL():
    file = open('raw_exp/JYFL2017.dat','r')
    lines = file.readlines()
    count = 0
    temp_d = {}
    for line in lines:
        if count<1: count=1; continue
        ss = line.split()
        Z,N = int(ss[0]),int(ss[1])
        me,meERR = float(ss[4]),float(ss[5])
        BE_ME = round(ME_BE(me,Z,N)/1000.0,6)
        meERR = round(meERR/1000.0,6)
        temp_d[(Z,N)] = (BE_ME, meERR)
    file.close()
    return temp_d

# returns dictionary with key (Z,N), value (BE, BE_ERR), energies are in MeV
def other_exp(name):
    file = open('raw_exp/'+name+'.dat','r')
    lines = file.readlines()
    count = 0
    temp_d = {}
    for line in lines:
        if count<1: count=1; continue
        ss = line.split()
        Z,N = int(ss[0]),int(ss[1])
        me,meERR = float(ss[2]),float(ss[3])
        BE_ME = round(ME_BE(me,Z,N)/1000.0,6)
        meERR = round(meERR/1000.0,6)
        temp_d[(Z,N)] = (BE_ME, meERR)
    file.close()
    return temp_d


# returns complete dictionary from Erik Olsen
# first layer key is experiment identifier, ie. '2003', 'J', etc.
# second layer key is (Z,N), value (BE, BE_ERR), energies are in MeV
def EO_data():
    # Create dictionary of EO's data.
    EO_dict = {}
    for k in ['2003','2016','J','T','R','O']:
            EO_dict[k] = {}
    for P1 in ['even','odd']:
        for P2 in ['even','odd']:
            f = open('Erik_BE/data_eoyc_2019_BE_'+P1+'_Z_'+P2+'_N.csv','r')
            l = f.readlines()
            count = 0
            for line in l:
                if count<1: count=1; continue
                ss = line.split(',')
                Z,N = int(ss[0]),int(ss[1])
                # EO doesn't have error for binding energy, store -1
                if ss[2] != '*':
                    EO_dict['2003'][(Z,N)] = (float(ss[2]), -1) #, float( ss[3]))
                if ss[4] != '*':
                    EO_dict['2016'][(Z,N)] = (float(ss[4]), -1) #, float( ss[5]))
                if ss[6] != '*':
                    EO_dict['J'][(Z,N)]    = (float(ss[6]), -1) #, float( ss[7])
                if ss[8] != '*':
                    EO_dict['T'][(Z,N)]    = (float(ss[8]), -1) #, float( ss[9])
                if ss[10] != '*':
                    EO_dict['R'][(Z,N)]    = (float(ss[10]), -1) #,float(ss[11])
                if ss[12] != '*':
                    EO_dict['O'][(Z,N)]    = (float(ss[12]), -1) #,float(ss[13])
    return EO_dict


def check(d1, d2):
    oStr1 = 'Z, N, experiment, who missed, value\n'
    oStr2 = 'Z, N, experiment, MC, EO, MC-EO, warning\n'
    for k in d1.keys():
        MC_nuc = set(d1[k].keys())
        EO_nuc = set(d2[k].keys())
        # Cross check missing nuclei
        M1 = MC_nuc - EO_nuc
        M2 = EO_nuc - MC_nuc
        common = MC_nuc & EO_nuc
        if M1:
            for Z,N in M1:
                if (Z,N) not in d2['J']:
                    oStr1 += str(Z)+','+str(N)+','+k+','+'Erik,'+str(d1[k][(Z,N)][0])+'\n'
                    #if d1[k][(Z,N)][0] > 0: print (Z,N,d1[k][(Z,N)][0])
        if M2:
            for Z,N in M2:
                oStr1 += str(Z)+','+str(N)+','+k+','+'Maxwell,'+str(d2[k][(Z,N)][0])+'\n'
        # Check common nuclei's difference
        for Z,N in common:
            dif = round(d1[k][(Z,N)][0] - d2[k][(Z,N)][0],6)
            oStr2 += str(Z)+','+str(N)+','+k+','+str(d1[k][(Z,N)][0])+','+\
            str(d2[k][(Z,N)][0])+','+str(dif)+','
            if abs(dif) >= 0.01:
                oStr2 += 'LARGE ERROR\n'
            else:
                oStr2 += '\n'
    return oStr1, oStr2

def hfb24_mass():
    f = open("raw_exp/hfb24.dat","r")
    l = f.readlines()
    f.close()
    temp_d = {}
    for line in l:
        ss = line.split()
        Z, A, me = int(ss[0]), int(ss[1]), float(ss[9])*1000.0
        N = A - Z
        temp_d[(Z,N)] = round(ME_BE(me,Z,N)/1000.0,2)
        #print (Z,N,temp_d[(Z,N)],me)
    return temp_d

def frdm_mass():
    f = open("raw_exp/FRDM2012.dat","r")
    l = f.readlines()
    f.close()
    temp_d = {}
    for line in l:
        ss = line.split()
        Z, N, be = int(ss[0]), int(ss[1]), float(ss[13])
        temp_d[(Z,N)] = be
        #print (Z,N,temp_d[(Z,N)],me)
    return temp_d

def qa_be(d,M_He):
    temp_d = {}
    for (Z,N) in d:
        if (Z-2,N-2) in d:
            temp_d[(Z,N)] = -round( d[(Z,N)] - d[(Z-2,N-2)] - M_He,2)
    return temp_d

def main(check_be=False,be_save=False,qa_save=False,check_qa=False):
    
    MC_dict, EO_dict = {}, EO_data()

    MC_dict['2003'] = AME_me2be(2003)
    MC_dict['2016'] = AME_me2be(2016)
    
    # TRIUMF, RIKEN, OTHERS, JYFLTRAP
    MC_dict['T'],MC_dict['R'] = other_exp('TRIUMF2018'),other_exp('RIKEN2018')
    MC_dict['O'] = other_exp('new_other')
    MC_dict['J'] = JYFL()
    
    M_He = MC_dict['2016'][(2,2)][0] # Mass of 4-He in MeV, from AME2016
    
    hfb24_be, frdm_be = hfb24_mass(), frdm_mass()
    hfb24_qa, frdm_qa = qa_be(hfb24_be,M_He), qa_be(frdm_be,M_He)
    
    if be_save:
        for P1 in ['even','odd']:
            for P2 in ['even','odd']:
                f = open('Erik_BE/data_eoyc_2019_BE_'+P1+'_Z_'+P2+'_N.csv','r')
                output = open('MC_BE/data_mcfixed_2019_BE_'+P1+'_Z_'+P2+'_N.csv','w')
                outStr = ''
                l = f.readlines()
                count = 0
                for line in l:
                    if count<1:
                        count=1; outStr = line
                        continue
                    ss = line.split(',')
                    Z,N = int(ss[0]),int(ss[1])
                    nuc = (Z,N)
                    nl_list = [str(Z),str(N)]
                    for k in ['2003','2016','J','T','R','O']:
                        if nuc in MC_dict[k]:
                            nl_list.extend([str(MC_dict[k][nuc][0]),str(MC_dict[k][nuc][1])])
                        else:
                            nl_list.extend(['*','*'])
                    nl_list.extend(ss[-10:-1])
                    # Erik's HFB-24 is wrong, updating as well
                    if nuc in hfb24_be:
                        nl_list.append(str(hfb24_be[(Z,N)]))
                    else:
                        nl_list.append('*')
                    nl = ','.join(nl_list)
                    outStr += nl+'\n'
                output.write(outStr)
                output.close()
    if qa_save:
        # Q_alpha = BE(Z,N) - BE(Z-2,N-2) - M(4He), BE is positive here
        MC_qa = {}
        MC_qa['2003'], MC_qa['2016'] = AME_qa(2003), AME_qa(2016)
        for k in ['J','T','R','O']:
            MC_qa[k] = {}
            for Z,N in MC_dict[k]:
                # Determine Q_alpha only for non-AME data, for AME data, use published data
                # Use most recent mass for (Z-2,N-2) by iterating the experiments in chrono. order
                # So that the mass get overwritten by the newest if there is any
                # If (Z-2,N-2) doesn't exist, make sure no qalpha value is saved
                # (Z,N) in non-AME set but (Z-2,N-2) in another set
                M_base = None
                for k_base in ['2003','2016','J','T','R','O']:
                    if (Z-2,N-2) in MC_dict[k_base]:
                        M_base = MC_dict[k_base][(Z-2,N-2)][0]
                        M_baseERR = MC_dict[k_base][(Z-2,N-2)][1]
                    # do not check experiment later than itself, it'll only be updated
                    # in the newer experiment's column if this happens
                    if k_base == k: break
                if M_base:
                    temp_qa = -round(MC_dict[k][(Z,N)][0] - M_base - M_He,6)
                    temp_qaERR = round(math.sqrt(MC_dict[k][(Z,N)][1]**2.0 + M_baseERR**2.0), 6)
                    if temp_qa>0: print ('daughter',k,Z,N,temp_qa,temp_qaERR)
                    MC_qa[k][(Z,N)] = (temp_qa,temp_qaERR)
                # (Z,N) in non-AME set but (Z+2,N+2) in another set
                M_parent = None
                for k_base in ['2003','2016','J','T','R','O']:
                    if (Z+2,N+2) in MC_dict[k_base]:
                        M_parent = MC_dict[k_base][(Z+2,N+2)][0]
                        M_parentERR = MC_dict[k_base][(Z+2,N+2)][1]
                    # do not check experiment later than itself, it'll only be updated
                    # in the newer experiment's column if this happens
                    if k_base == k: break
                if M_parent:
                    temp_qa = -round(M_parent - MC_dict[k][(Z,N)][0] - M_He,6)
                    temp_qaERR = round(math.sqrt(MC_dict[k][(Z,N)][1]**2.0 + M_parentERR**2.0), 6)
                    if temp_qa>0: print ('parent',k,Z+2,N+2,temp_qa,temp_qaERR)
                    MC_qa[k][(Z+2,N+2)] = (temp_qa,temp_qaERR)

        EO_qa = {}
        for k in ['2003','2016','J','T','R','O']:
            EO_qa[k] = {}

        for P1 in ['even','odd']:
            for P2 in ['even','odd']:
                f = open('Erik_Qalpha/data_eoyc_2019_Qalpha_'+P1+'_Z_'+P2+'_N.csv','r')
                output = open('MC_Qalpha/data_mcfixed_2019_Qalpha_'+P1+'_Z_'+P2+'_N.csv','w')
                outStr = ''
                l = f.readlines()
                count = 0
                for line in l:
                    if count<1:
                        count=1; outStr = line
                        continue
                    ss = line.split(',')
                    Z,N = int(ss[0]),int(ss[1])
                    nuc = (Z,N)
                    nl_list = [str(Z),str(N)]
                    for ind,k in enumerate(['2003','2016','J','T','R','O']):
                        temp_qa, temp_qaERR = ss[2*(ind+1)], ss[2*(ind+1)+1]
                        if temp_qa != '*':
                            EO_qa[k][(Z,N)] = (float(temp_qa),float(temp_qaERR))
                    for k in ['2003','2016','J','T','R','O']:
                        if nuc in MC_qa[k]:
                            nl_list.extend([str(MC_qa[k][nuc][0]),str(MC_qa[k][nuc][1])])
                        else:
                            nl_list.extend(['*','*'])
                    nl_list.extend(ss[-10:-2])
                    if nuc in frdm_qa:
                        nl_list.append(str(frdm_qa[(Z,N)]))
                    else:
                        nl_list.append('*')
                    if nuc in hfb24_qa:
                        nl_list.append(str(hfb24_qa[(Z,N)]))
                    else:
                        nl_list.append('*')
                    nl = ','.join(nl_list)
                    outStr += nl+'\n'
                output.write(outStr)

    if check_be:
        o1 = open('EO_MC_be_missing.csv','w')
        o2 = open('EO_MC_be_diff.csv','w')
        oStr1, oStr2 = check(MC_dict, EO_dict)
        o1.write(oStr1); o2.write(oStr2)
        o1.close(); o2.close()

    if check_qa:
        o1 = open('EO_MC_qa_missing.csv','w')
        o2 = open('EO_MC_qa_diff.csv','w')
        oStr1, oStr2 = check(MC_qa, EO_qa)
        o1.write(oStr1); o2.write(oStr2)
        o1.close(); o2.close()

main(qa_save=True,check_qa=True)



        








