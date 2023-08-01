# IPython log file
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
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
def f_2sq(s, prof):
    l = len(prof)
    ls = int(l/s)
    v_ls = []
    for v in range(ls):
        i0 = int(v*s)
        iff = int((v+1)*s)
        cur_ls = prof[i0:iff]
        x = list(range(len(cur_ls)))
        y_curl = np.polyfit(x, cur_ls, 3) # give 4 values
        y1, y2, y3, y4 = y_curl
        cur_v = 0
        for i in range(len(cur_ls)):
            ii = (cur_ls[i] - y1*(i**3) - y2*(i**2) - y3*i - y4)**2
            cur_v = cur_v + ii
        v_ls.append(cur_v/s)
    return v_ls

def f_q(f2, s, q):
    ls = len(f2)
    allf = 0
    for f in f2:
        ff = f**(q/2)
        allf = allf + ff
    result = (1/ls*allf)**(1/q)
    return result

def hs_con(fq_ls, s_ls):
    fq_ln = []
    s_ln = []
    for i in range(len(s_ls)):
        fq_ln.append(np.log(fq_ls[i]))
        s_ln.append(np.log(s_ls[i]))
    gr, con = np.polyfit(s_ln, fq_ln, 1)
    return gr
def tau_ls(hs, q):
    return hs*q - 1

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

def extract_file(file, si, fi, ex_id):
    # 152, 1152, 2
    sf = fi - si
    mylines = []
    with open(file, 'rt') as myfile:
        for myline in myfile:
            mylines.append(myline)
    x = range(0, sf)
    real_num = []
    for i in range(si, fi):
        iv = mylines[i].split()[ex_id]
        ii = round(eval(iv), 4)
        real_num.append(ii)
    return real_num
def profile_ls(real_num):
    avg = sum(real_num)/len(real_num)
    profile = []
    cur_sum = 0
    for i in real_num:
        iele = i - avg + cur_sum
        profile.append(iele)
        cur_sum = iele
    return profile

def give_f_z(profile, q_ls, s_ls):
    all_tau = []
    for q in q_ls:
        cnt_fq = []
        for s in s_ls:
            f2 = f_2sq(s, profile)
            fq = f_q(f2, s, q)
            cnt_fq.append(fq)
        hs = hs_con(cnt_fq, s_ls)
        tau = tau_ls(hs, q)
        all_tau.append(tau)
    al_ls = alpha_fd(q_ls, all_tau)
    f_ls = []
    for a in range(len(al_ls)):
        aa = a + 1
        ff = al_ls[a]*q_ls[aa] - all_tau[aa]
        f_ls.append(ff)
        #print(ff, al_ls[a], q_ls[aa])
    return f_ls, al_ls
    
#q_ls = list(np.arange(-2, 1, 0.3))
s_ls = list(np.arange(5, 19, 2))
k = 't02.lammpstrj'
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
sidx = nline + 6
eidx = sidx + natom
        
        
q_ls = list(np.arange(-2, 2, 0.3))
    
def partition(unitlen, q):
    brange = [0, 20]
    #unitlen = 4
    tns = int((abs(brange[1] - brange[0])/unitlen)**3)
    ttlist = np.zeros(tns)
    for j in range(sidx, eidx):
        coor = raw_line(mylines[j])[2]
        boxidx = int(ide_box(brange, unitlen, coor))
        if boxidx < len(ttlist):
            ttlist[boxidx] = ttlist[boxidx] + 1
        #else:
        #    ttlist[boxidx] = ttlist[boxidx] + 1
    pp = (ttlist/(1*10**(-8)+sum(ttlist)))**q
    return pp
   
#print(q_ls)
singleq = [] # this is tau list
q_lsr = []
for q in q_ls:
    allcp = []
    for s in s_ls:
        nn = np.log(sum(partition(s, q)))/np.log(s/20)
        #print(nn)
        if np.round(nn, 2) != 0 and nn != math.inf and nn != -math.inf:
            allcp.append(nn)
    #print(allcp)
    if len(allcp) == 0:
        q_lsr.append(q)
    else:
        singleq.append(sum(allcp)/len(allcp))
    #print('------------')
for qr in q_lsr:
    q_ls.remove(qr)
#print('the new q')
#print(q_ls)
#print(singleq)
alpha = alpha_fd(q_ls, singleq)
f_ls = f_fd(alpha, singleq, q_ls)
#print(f_ls)

print(q_ls)
print(singleq)
print(alpha)
print(f_ls)

plt.plot(alpha, f_ls)
plt.xlabel('Singularity Strength')
plt.ylabel('Singularity Spectrum')
plt.title('Multifractal Spectrum')
#plt.legend()
plt.show()

