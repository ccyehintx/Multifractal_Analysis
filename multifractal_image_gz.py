# IPython log file
import gzip
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import json
def raw_line(line):
    info = []
    info.append(eval(line.split()[0]))
    info.append(eval(line.split()[2]))
    coor = np.array([eval(i) for i in line.split()[3:6]])
    info.append(coor)
    return info

def sort_info_func(mylines, sidx, fidx):
    sort_info = []
    for i in range(sidx+9, fidx+1):
        sort_info.append(raw_line(mylines[i]))
    return sort_info

def colloid_ls(sort_info, colloid_idx):
    c_ls = []
    for i in range(len(sort_info)):
        if sort_info[i][1] == colloid_idx:
            c_ls.append(sort_info[i])
    return c_ls

def ide_box(brange, unitlen, coor):
    nsubs = (abs(brange[1] - brange[0])/unitlen)
    tns = nsubs**3
    x,y,z = coor - np.array([0.001, 0.001, 0.001])
    zbox = (z - brange[0])//unitlen
    ybox = (y - brange[0])//unitlen
    xbox = (x - brange[0])//unitlen
    idx = xbox + ybox*nsubs + zbox*nsubs**2
    return idx

def alpha_fd(x, y):
    # apply finite difference method
    # y is tau list and x is q list
    alpha_ls = []
    for i in range(len(x)-1):
        x0, x1 = x[i], x[i+1]
        y0, y1 = y[i], y[i+1]
        alpha = (y1 - y0)/(x1 - x0)
        alpha_ls.append(alpha)
        # grab all q except the first one
    return alpha_ls

def f_fd(alpha, tau, q):
    a_ls = alpha
    t_ls = tau
    qq_ls = q
    del t_ls[0]
    del qq_ls[0]
    f_ls = []
    for i in range(len(a_ls)):
        f = a_ls[i]*qq_ls[i] - t_ls[i]
        f_ls.append(f)
    return f_ls

def partition(brange, unitlen, q, coor_ls):
    tns = int((abs(brange[1] - brange[0])/unitlen)**3)
    ttlist = np.zeros(tns)
    for j in coor_ls:
        coor = j[2]
        boxidx = int(ide_box(brange, unitlen, coor))
        if boxidx < len(ttlist):
            ttlist[boxidx] = ttlist[boxidx] + 1
    pp = (ttlist/(1*10**(-8)+sum(ttlist)))**q
    return pp

f = open('input_gz.json')
data = json.load(f)
gzfile = data['file_lists']['filename']
#filels = data['file_lists']['filename']
lp_co = eval(data['specifications']['line_plus'])
exsnap = data['specifications']['extractsnap']
exincr = eval(data['specifications']['extract_incre'])
brange = data['specifications']['box_range']
sls0 = eval(data['specifications']['box_slice_0'])
slsf = eval(data['specifications']['box_slice_f'])
slsd = eval(data['specifications']['box_slice_d'])
qls0 = eval(data['specifications']['q_list_0'])
qlsf = eval(data['specifications']['q_list_f'])
qlsd = eval(data['specifications']['q_list_d'])
q_ls = list(np.arange(qls0, qlsf, qlsd))
s_ls = list(np.arange(sls0, slsf, slsd))

mylines = []
with gzip.open(gzfile, mode ='r')as file:
    for ll in file:
        mylines.append(ll)

for nstep in range(exsnap[0], exsnap[1], exincr):
    qc_ls = q_ls.copy()
    sidx = 31009*nstep
    eidx = sidx + 31008
    sort_info = sort_info_func(mylines, sidx, eidx)
    c_ls = colloid_ls(sort_info, 1) # colloid idx
    # Calculation of multifractal begins
    singleq = [] # this is tau list
    q_lsr = []
    for q in qc_ls:
        allcp = []
        for s in s_ls:
            nn = np.log(sum(partition(brange, s, q, c_ls)))/np.log(s/(brange[1] - brange[0]))
            if np.round(nn, 2) != 0 and nn != math.inf and nn != -math.inf:
                allcp.append(nn)
        if len(allcp) == 0:
            q_lsr.append(q)
        else:
            singleq.append(sum(allcp)/len(allcp))
    for qr in q_lsr:
        qc_ls.remove(qr)
    alpha = alpha_fd(qc_ls, singleq)
    print(nstep, (max(alpha) - min(alpha)))
    f_ls = f_fd(alpha, singleq, qc_ls)
    ffname = nstep
    plt.plot(alpha, f_ls, label=ffname)

plt.xlabel('Singularity Strength')
plt.ylabel('Singularity Spectrum')
plt.title('Multifractal Spectrum')
plt.legend()
plt.show()

