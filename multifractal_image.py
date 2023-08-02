# IPython log file
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import json
def raw_line(line):
    info = []
    info.append(eval(line.split(' ')[0]))
    info.append(eval(line.split(' ')[1]))
    coor = np.array([eval(i) for i in line.split(' ')[2:5]])
    info.append(coor)
    return info


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

def partition(brange, unitlen, q):
    tns = int((abs(brange[1] - brange[0])/unitlen)**3)
    ttlist = np.zeros(tns)
    for j in range(sidx, eidx):
        coor = raw_line(mylines[j])[2]
        boxidx = int(ide_box(brange, unitlen, coor))
        if boxidx < len(ttlist):
            ttlist[boxidx] = ttlist[boxidx] + 1
    pp = (ttlist/(1*10**(-8)+sum(ttlist)))**q
    return pp

f = open('input.json')
data = json.load(f)
filels = data['file_lists']['filename']
lp_co = eval(data['specifications']['line_plus'])
brange = data['specifications']['box_range']
sls0 = eval(data['specifications']['box_slice_0'])
slsf = eval(data['specifications']['box_slice_f'])
slsd = eval(data['specifications']['box_slice_d'])
qls0 = eval(data['specifications']['q_list_0'])
qlsf = eval(data['specifications']['q_list_f'])
qlsd = eval(data['specifications']['q_list_d'])
q_ls = list(np.arange(qls0, qlsf, qlsd))
s_ls = list(np.arange(sls0, slsf, slsd))

for k in filels:
    filename = k
    mylines = []
    with open(filename, mode ='r')as file:
        for ll in file:
            mylines.append(ll)
    cmd = "grep -n 'ITEM: NUMBER OF ATOMS' {}".format(filename)
    command = subprocess.check_output(cmd, shell=True)
    cc = command.decode("utf-8").split("\n")
    i = cc[50]
    nline = eval(i.split(':')[0])
    natom = eval(mylines[nline])
    sidx = nline + lp_co
    eidx = sidx + natom
    
    singleq = [] # this is tau list
    q_lsr = []
    for q in q_ls:
        allcp = []
        for s in s_ls:
            nn = np.log(sum(partition(brange, s, q)))/np.log(s/(brange[1] - brange[0]))
            if np.round(nn, 2) != 0 and nn != math.inf and nn != -math.inf:
                allcp.append(nn)
        if len(allcp) == 0:
            q_lsr.append(q)
        else:
            singleq.append(sum(allcp)/len(allcp))
    for qr in q_lsr:
        q_ls.remove(qr)
    alpha = alpha_fd(q_ls, singleq)
    f_ls = f_fd(alpha, singleq, q_ls)
    ffname = k.split('.')[0]
    plt.plot(alpha, f_ls, label=ffname)

plt.xlabel('Singularity Strength')
plt.ylabel('Singularity Spectrum')
plt.title('Multifractal Spectrum')
plt.legend()
plt.show()

