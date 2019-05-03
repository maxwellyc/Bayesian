import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import math
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_file( file_loc ):
    f1 = open(file_loc)
    l1 = f1.readlines()
    S1n = {}; S2n = {}
    # S1/2n keys: (proton, neutron, models, S?n value type)
    # Models ordering:
    #           1       2     3     4     5     6      7       8       9
    #Leo's: FRDM2012, HFB24, SKM*, SKP, SLY4, SVMIN, UNEDF0, UNEDF1, UNEDF2
    
    #        1       2       3      4     5    6      7       8        9
    #def # UNEDF0, UNEDF1, UNEDF2, SLY4, SKP, SKM*, SVMIN, FRDM2012, HFB24
    
    # S1/2n value types (0~13):
    #  0     1      2      3        4       5       6       7       8       9      10      11       12      13
    # raw    gp    sd    0.005    0.025    0.05    0.16    0.33    0.5    0.67    0.84    0.95    0.975    0.995
    # EXCEPTION: SKM* S2n value is missing 0.95, 0.975, 0.995.
    # reorder_20** maps current model ordering from input to default ordering of models
    reorder_2003 = [8,9,6,5,4,7,1,2,3]
    reorder_2016 = [8,9,6,5,4,7,1,2,3]
    reorder_2018 = [8,9,6,5,4,7,1,2,3]
    
    if "2003" in file_loc:
        for line in l1:
            # remove last char "\n" for change line
            line = line[0:-1]
            #if the predictions_*.dat is produced by copy and pasting macOS numbers, use "    " (4 space) separator
            ss = line.split("    ")
            #if the predictions_*.dat is produced by copy and pasting windows excel, use "\t" separator
            #ss = line.split("\t")
            ss_len = len(ss)
            Z = int(ss[0])
            N = int(ss[1])
            A = int(ss[2])
            if N%2 == 1:
                if len(ss[3]) > 0:  S2n[(Z,N,0,0)] = round(float(ss[3]),6)
                for i in range(1,10):
                    ii = reorder_2003[i-1]
                    for j in range(0,14):
                        if (5+(i-1)*14+j) < ss_len:
                            v_str = ss[5+(i-1)*14+j]
                            if len(v_str) != 0:
                                S1n[(Z,N,ii,j)] = round(float(v_str),6)
            elif N%2 == 0:
                if len(ss[4]) > 0:  S1n[(Z,N,0,0)] = round(float(ss[4]),6)
                for i in range(1,10):
                    ii = reorder_2003[i-1]
                    for j in range(0,14):
                        if (131+(i-1)*14+j) < ss_len:
                            v_str = ss[131+(i-1)*14+j]
                            if len(v_str) != 0:
                                S2n[(Z,N,ii,j)] = round(float(v_str),6)
        return S1n, S2n
    
    elif "2016" in file_loc:
        for line in l1:
            # remove last char "\n" for change line
            line = line[0:-1]
            #if the predictions_*.dat is produced by copy and pasting macOS numbers, use "    " (4 space) separator
            ss = line.split("    ")
            #if the predictions_*.dat is produced by copy and pasting windows excel, use "\t" separator
            #ss = line.split("\t")
            ss_len = len(ss)
            Z = int(ss[0])
            N = int(ss[1])
            A = int(ss[2])
            if N%2 == 1:
                if len(ss[3]) > 0:  S1n[(Z,N,0,0)] = round(float(ss[3]),6)
                for i in range(1,10):
                    ii = reorder_2016[i-1]
                    for j in range(0,14):
                        if (5+(i-1)*14+j) < ss_len:
                            v_str = ss[5+(i-1)*14+j]
                            if len(v_str) != 0:
                                S1n[(Z,N,ii,j)] = round(float(v_str),6)
            elif N%2 == 0:
                if len(ss[4]) > 0:  S2n[(Z,N,0,0)] = round(float(ss[4]),6)
                for i in range(1,10):
                    ii = reorder_2016[i-1]
                    for j in range(0,14):
                        if (131+(i-1)*14+j) < ss_len:
                            v_str = ss[131+(i-1)*14+j]
                            if len(v_str) != 0:
                                S2n[(Z,N,ii,j)] = round(float(v_str),6)
        return S1n, S2n

    if "2018" in file_loc:
        for line in l1:
            # remove last char "\n" for change line
            line = line[0:-1]
            #if the predictions_*.dat is produced by copy and pasting macOS numbers, use "    " (4 space) separator
            ss = line.split("    ")
            #if the predictions_*.dat is produced by copy and pasting windows excel, use "\t" separator
            #ss = line.split("\t")
            ss_len = len(ss)
            Z = int(ss[0])
            N = int(ss[1])
            A = int(ss[2])
#            print (ss_len)
#            print (ss)
            if N%2 == 1:
                if len(ss[3]) > 0:  S1n[(Z,N,0,0)] = round(float(ss[3]),6)
                for i in range(1,10):
                    ii = reorder_2018[i-1]
                    for j in range(0,14):
                        if (5+(i-1)*14+j) < ss_len:
                            v_str = ss[5+(i-1)*14+j]
                            if len(v_str) != 0:
                                S1n[(Z,N,ii,j)] = round(float(v_str),6)
            elif N%2 == 0:
                if len(ss[4]) > 0:  S2n[(Z,N,0,0)] = round(float(ss[4]),6)
                for i in range(1,10):
                    ii = reorder_2018[i-1]
                    for j in range(0,14):
                        if (131+(i-1)*14+j) < ss_len:
                            v_str = ss[131+(i-1)*14+j]
                            if len(v_str) != 0:
                                S2n[(Z,N,ii,j)] = round(float(v_str),6)
        return S1n, S2n

# proton#: plot_z and neutron#: plot_n
def plot_sep_model():
    plot_z = 20; Tg = 'Ca'
    predictions_03 = read_file("predictions_2003_new.dat")
    S1n_03, S2n_03 = predictions_03
    predictions_16 = read_file("predictions_2018_new_noRIKEN.dat")
    S1n_16, S2n_16 = predictions_16
    print (S1n_16[(20,35,9,0)])
    print (S2n_16[(20,36,9,0)])
    print (S1n_16[(20,37,9,0)])
    model_name = ["exp","UNEDF0", "UNEDF1", "UNEDF2", "SLy4", "SkP", "SkM*", "SV-min", "FRDM-2012", "HFB-24" ]
    plt.clf()
    dis = 0.2 # x direction distance between same functional, different predictions
#    ax = []
#    ax.append(plt.subplot(212))
#    ax.append(plt.subplot(221))
#    ax.append(plt.subplot(222))
    rows = 13; cols = 2
    fig, ax = plt.subplots(rows, cols, sharex = True, sharey = True)#, figsize = (30, 45))
    gs = gridspec.GridSpec(rows, cols, wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.17, right=0.845)
    ax[0,0] = plt.subplot(gs[0:7,0])  #55Ca
    ax[1,0] = plt.subplot(gs[:,1])   #56Ca
    ax[1,1] = plt.subplot(gs[7:13,0])  #57Ca
    for plot_n in range(35,38):
        if plot_n == 35: p = 0; q = 0
        elif plot_n == 36: p = 1; q = 0
        elif plot_n == 37: p = 1; q = 1
        exp_ind = plot_n - 35 # ordering of below exp measurements, to access *_meas array and its *_sd
        #55Ca, 56Ca, 57Ca, S1n for odd N, S2n for even N
        new_meas = [1.560733, 4.492051, 1.931318]
        new_meas_sd = [0.167171, 0.254649, 1.021078]
        fig_name = r"$^{55}$Ca"
        ax[p,q].set_title(fig_name)
        x_exp_new = []; y_exp_new = []; y_sd_exp_new = []
        x_exp_old = []; y_exp_old = []; y_sd_exp_old = []
        x_ticks = []; y_raw = []; x_raw = []
        y_03 = []; y_sd_03 = []; x_03 = [];
        y_16 = []; y_sd_16 = []; x_16 = [];
        # xticks in order of: UNEDF0~2, SLy4, SkP, SkM*, SV-min, FRDM-2012, HFB-24
        for i in range(0,10):
            x_ticks.append(model_name[i])
        # Add new measurement 2018
        x_exp_new.append(0);
        y_exp_new.append(new_meas[exp_ind])
        y_sd_exp_new.append(new_meas_sd[exp_ind])
        # Add older measurement
        if plot_z == 22:
            x_exp_old.append(0);
            y_exp_old.append(old_meas[exp_ind])
            y_sd_exp_old.append(old_meas_sd[exp_ind])
        if not plot_n%2:
            count = 0
            for i in range(1,10):
                count = count + 1
                if (plot_z,plot_n,i,1) in S2n_03:
                    #03 prediction at x-dis
                    x_03.append(count-dis)
                    y_03.append(S2n_03[(plot_z,plot_n,i,1)])
                    y_sd_03.append(S2n_03[(plot_z,plot_n,i,2)])
                    #raw at x
                    x_raw.append(count)
                    y_raw.append(S2n_03[(plot_z,plot_n,i,0)])
                    #16/18 prediction at x+dis
                    x_16.append(count+dis)
                    y_16.append(S2n_16[(plot_z,plot_n,i,1)])
                    y_sd_16.append(S2n_16[(plot_z,plot_n,i,2)])
        elif plot_n%2:
            count = 0
            for i in range(1,10):
                count = count + 1
                if (plot_z,plot_n,i,1) in S1n_03:
                    #03 prediction at x-dis
                    x_03.append(count-dis)
                    y_03.append(S1n_03[(plot_z,plot_n,i,1)])
                    y_sd_03.append(S1n_03[(plot_z,plot_n,i,2)])
                    #raw at x
                    x_raw.append(count)
                    y_raw.append(S1n_03[(plot_z,plot_n,i,0)])
                    #16 prediction at x+dis
                    x_16.append(count+dis)
                    y_16.append(S1n_16[(plot_z,plot_n,i,1)])
                    y_sd_16.append(S1n_16[(plot_z,plot_n,i,2)])
        y_03 = np.array(y_03)
        y_sd_03 = np.array(y_sd_03)
        y_16 = np.array(y_16)
        y_sd_16 = np.array(y_sd_16)
        y_raw = np.array(y_raw)
        ax[p,q].set_xticks(np.arange(10), x_ticks)
        ax[p,q].errorbar(x_exp_new, y_exp_new, yerr=y_sd_exp_new, color = 'k', ecolor = 'k', fmt='o',capsize=5, label = 'exp_2018')
        if plot_n == 35:
            ax[p,q].set_ylim([0.5,4])
            ax[p,q].set_yticks(np.arange(0.5,4.1,0.5))
            #AME2016 extrapolated value:
            ax[p,q].errorbar([0.2], [1.261], yerr=[0.304], color = 'w', ecolor = 'k',mec='k', fmt='s',capsize=5, \
                 label = 'AME2016 extrapolation')
        elif plot_n == 36:
            ax[p,q].set_ylim([1,7.5])
            ax[p,q].set_yticks(np.arange(1,8,0.5))
            #AME2016 extrapolated value:
            ax[p,q].errorbar([0.2], [4.88], yerr=[0.403], color = 'w', ecolor = 'k',mec='k', fmt='s',capsize=5, \
                 label = 'AME2016 extrapolation')
        elif plot_n == 37:
            ax[p,q].set_ylim([0,3])
            ax[p,q].set_yticks(np.arange(0,3.1,0.5))
            #AME2016 extrapolated value:
            ax[p,q].errorbar([0.2], [1.048], yerr=[0.565], color = 'w', ecolor = 'k',mec='k', fmt='s',capsize=5, \
                 label = 'AME2016 extrapolation')
    #    if plot_z == 22:
    #        plt.errorbar(x_exp_old, y_exp_old, yerr=y_sd_exp_old, color = 'b', ecolor = 'b', fmt='x',capsize=5, label = 'exp_2016')
    #
        exp_ll = new_meas[exp_ind] - new_meas_sd[exp_ind]
        exp_hl = new_meas[exp_ind] + new_meas_sd[exp_ind]
        ax[p,q].fill([0,0,9.2,9.2],[exp_ll-0.005,exp_hl+0.005,exp_hl+0.005,exp_ll-0.005], \
            color = 'gray', alpha=0.2, facecolor='none', lw=0)
    #    plt.plot([0,9.1],[exp_ll, exp_ll], 'k--', lw = 1)
    #    plt.plot([0,9.1],[exp_hl, exp_hl], 'k--', lw = 1)
        ax[p,q].plot([0,9.2],[new_meas[exp_ind], new_meas[exp_ind]], 'k--', lw = 0.5)
        ax[p,q].errorbar(x_03,y_03,yerr=y_sd_03, fmt='^', color = 'w', mec='r',ecolor = 'r', capsize=5, label = r'gp_03+$\sigma$')
        ax[p,q].errorbar(x_16,y_16,yerr=y_sd_16, fmt='v', color = 'b', ecolor = 'b', capsize=5, label = r'gp_18+$\sigma$')
        #if exp_ind == 0:
        #    plt.errorbar(sbnn_x, sbnn_y, yerr = sbnn_sd, color = 'c', ecolor = 'c', capsize=5, label = 'model+sbnn', fmt = '*')
        # plot raw model
        ax[p,q].scatter(x_raw, y_raw, c = '#00FF00', label = 'raw_model', marker = 'd', s = 50, edgecolor='k')
        if plot_n == 36: plt.legend(); ax[p,q].yaxis.tick_right()
    fig = plt.gcf()
    fig.set_size_inches(15,10)
    fig.savefig('55_57Ca_03_18.pdf', dpi=100,format = "pdf")


    #plt.show()

def plot_delta_2n():
    predictions = read_file("predictions_data.dat")
    S1n, S2n = predictions
    #model_name = ["exp","UNEDF1", "FRDM-2012", "UNEDF2", "SLY4", "SVMIN", "UNEDF0", "SKP", "HFB-24", "SKM*"]
    model_name = ["exp","SkP", "HFB-24", "UNEDF0", "SLy4", "SkM*", "FRDM-2012", "UNEDF1", "UNEDF2", "SV-min"]
    delta_2n = {}; delta_2n_sd = {}
    x_d2n = []; y_d2n = []
    Z = 22
    x_d2n.append([])
    y_d2n.append([])
    for N in range(12,60,2):
        if (Z,N,0,0) in S2n and (Z,N+2,0,0) in S2n:
            delta_2n[(Z,N,0)] = round(S2n[(Z,N,0,0)] - S2n[(Z,N+2,0,0)],6)
            x_d2n[0].append(N)
            y_d2n[0].append(delta_2n[(Z,N,0)])
    for i in range(1,10):
        x_d2n.append([])
        y_d2n.append([])
        for N in range(12,60,2):
            if (Z,N,i,1) in S2n and (Z,N+2,i,1) in S2n:
                delta_2n[(Z,N,i)] = round(S2n[(Z,N,i,1)] - S2n[(Z,N+2,i,1)],6)
                x_d2n[i].append(N)
                y_d2n[i].append(delta_2n[(Z,N,i)])
                #delta_2n_sd[(Z,N,i)] = round(math.sqrt( (S2n[(Z,N,i,2)])**2 + (S2n[(Z,N+2,i,2)])**2 ),6)
#    print (x_d2n)
#    print (y_d2n)
    for i in range(0,10):
        if i == 0:
            mk = 'x'; ms = 5; colors = 'r'; linew = 1; lines = '--'
            plt.plot(np.array(x_d2n[i]), np.array(y_d2n[i]) , marker = mk, markersize = ms, c = colors, ls = lines, lw = linew, label = model_name[i])
        else:
            mk = '.'; ms = 1; linew = 1; lines = '-'
            plt.plot(np.array(x_d2n[i]), np.array(y_d2n[i]) , marker = mk, markersize = ms, ls = lines, lw = linew, label = model_name[i])
    plt.title(r"Ti $\Delta_{2n}$ vs N")
    plt.xlabel("N")
    plt.ylabel(r"$\Delta_{2n}$ (MeV)")
    plt.xticks(np.arange(10,60,5))
    plt.legend()
    plt.show()

#plot_sep_model(20,35)
plot_sep_model()
#plot_sep_model(20,37)
#plot_delta_2n()
