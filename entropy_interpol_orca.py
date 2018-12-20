#! /usr/bin/python

import sys, os
from optparse import OptionParser
from math import exp, log, expm1, log1p, pi, pow, sqrt
import matplotlib.pyplot as plt
import moment_inert
import molpro

homedir = "/home/stvogt/"

h = 6.62606957e-34
R = 8.3144621
Rkcal = 1.987191683e-3
T = 298.15
#T = 100.0
k = 1.3806488e-23
c = 29979245800
alpha = 4
w0 = 100
amuA2mkg2 = 1.6605e-47

parser = OptionParser() 
parser.add_option("-m", "--method", dest="method", help="The level of theory", default = "B3PW91") #, metavar="FILE")
parser.add_option("-b", "--basis", dest="basis", help="The basis set level", default="def2-tzvpp")
parser.add_option("-o", "--output", dest="output", help="The name of the output file", default="output.log")
parser.add_option("-u", "--units", dest="units", help="The conversion units", default="kcal")
parser.add_option("-r", "--reference", action="store_true", dest="ref", help=" The reference with respect reactant complex?")

(options, args) = parser.parse_args()

method = options.method
basis = options.basis
out = options.output
units = options.units
ref = options.ref

if units == 'kcal':
    conv = 627.509469
    convS = "kcal/mol"
if units == 'kj':
    conv = 2625.5
    convS = "kJ/mol"

def extr_energy(output):
    energy_dict = {}
    L_info = []
    try:
        lines = open(output,'r').readlines()
    except IOError:
        print "\n\nError!!   No output found in:  "+os.getcwd()
        sys.exit(1)
    for line in lines:
        #if "Electronic energy      " in line:
        #    en_el = line.split()[3]
        #    energy_dict["E_el"] = en_el
        #if "Final Gibbs free enthalpy      " in line:
        #    en_G = line.split()[5]
        #    energy_dict["en_G"] = en_G
        #if "G-E(el)  " in line:
        #    G_correct = line.split()[2]
        #    energy_dict["G_correct"] = G_correct
        #if "Total thermal correction  " in line:
        #    H_correct = line.split()[3]
        #    energy_dict["Therm_correct"] = H_correct
        if "Total thermal correction  " in line:
            Therm_correct = line.split()[3]
            energy_dict["Therm_correct"] = Therm_correct
        if "Thermal Enthalpy correction  " in line:
            ThermH_correct = line.split()[4]
            energy_dict["ThermH_correct"] = ThermH_correct
        if "Non-thermal (ZPE) correction " in line:
            ZPVE_correct = line.split()[3]
            energy_dict["ZPVE_correct"] = ZPVE_correct
        #if "Total entropy correction " in line:
        #    S_correct = line.split()[4]
        #    energy_dict["S_correct"] = S_correct
        #    #print "ENTROPY ======== "+S_correct
        #if "Rotational entropy   " in line:
        #    Srot_correct = line.split()[3]
        #    energy_dict["Srot_correct"] = Srot_correct
        if "Vibrational entropy   " in line:
            Svib_correct = line.split()[3]
            energy_dict["Svib_correct"] = Svib_correct
        #if "Translational entropy   " in line:
        #    Stran_correct = line.split()[3]
        #    energy_dict["Stran_correct"] = Stran_correct
    return energy_dict

def extr_entropy(output):
    entropy_list = []
    try:
        lines = open(output,'r').readlines()
    except IOError:
        print "\n\nError!!   No output found in:  "+os.getcwd()
        sys.exit(1)
    for lineNum in range(0,len(lines)):
        line = lines[lineNum]
        if "VIBRATIONAL FREQUENCIES" in line:
            for lineNum1 in range(lineNum,len(lines)):
                if "cm**-1" in lines[lineNum1]:
                    freq = lines[lineNum1].split()[1]
                    if float(freq) > 0.00:
                        entropy_list.append(float(freq))
                        #print lines[lineNum1].split()[1]
                if "NORMAL MODES" in lines[lineNum1]:
                    break
    return entropy_list



def Sv(v):
    w = float(v)*c
    hwkt=(h*w)/(k*T)
    tmp1 = hwkt/(expm1(hwkt))
    tmp2 =log1p(-exp(-hwkt))
    tmp=tmp1-tmp2
    S_v = Rkcal*(tmp1 -tmp2)
    return S_v

def mu(v):
    w = float(v)*c
    I = h/(8*pow(pi,2)*w)
    #print "mu = "+str(I)
    return I


def Bav(output):
    coords = moment_inert.read_coords_orca(output)
    COM_coords = moment_inert.center_of_mass(coords)
    I_list = moment_inert.moment_inertia(COM_coords)
    #print I_list
    sum_I = 0.0
    for I in I_list:
        sum_I = sum_I + I
    Bav = (sum_I/3)*amuA2mkg2
    #print Bav
    return Bav

    

#def Bav(freq_list):
#    count = 0
#    muav = 0.0
#    for freq in freq_list:
#        if float(freq) > 0.00:
#            muav = muav + mu(freq)
#            count = count + 1
#    B_av = muav/count
#    return B_av
    
def eff_mu(mu,Bav):
    I_eff = (mu*Bav)/(mu+Bav)
    #print "mu_prime: "+str(I_eff)
    return I_eff

def Sr(eff_mu):
    S_r = Rkcal*(0.5+log(sqrt((8*pow(pi,3)*eff_mu*k*T)/pow(h,2))))
    #print sqrt((8*pow(pi,3)*eff_mu*k*T)/pow(h,2))
    return S_r

def w(v):
    return 1/(1+pow((w0/float(v)), alpha))

runDir = os.getcwd()
all_energies = {}
r_en = {}
direc_list = []
method_list = [] 
results1 = {}
results2 = {}
cat_en = {}
cat_mon_en = {}
react_en = {}

print "\n\nMethod:  "+method
print "Basis:   "+basis 
print '{:20s}  {:^25s}  {:^38s}'.format('Structure(Method)',  'Relative Energies (in '+convS+')' , 'Total Energies')
dirs = os.listdir(runDir)
dir_dict = {}
dir_list = []

for d in dirs:
    if "done" in d:
        dir_dict[d.split('_')[0]] = d
        dir_list.append(d.split('_')[0])

if ref:
    print sorted(dir_list, key=molpro.struct_key_r_first)
    iter_struct = sorted(dir_list, key=molpro.struct_key_r_first)
else:
    print sorted(dir_list, key=molpro.struct_key_mon_first)
    iter_struct = sorted(dir_list, key=molpro.struct_key_mon_first)

for d in iter_struct:
    #d = dir_dict[d_tmp]
    try:
        os.chdir(dir_dict[d]+'/'+method+'/'+basis+'/freq')
        print '\n****************************************************************************\nProcessing structure: '+dir_dict[d]
    except OSError:
        print "\n\nThis directory does not exsit ----> " +os.getcwd()+'/'+d+'/'+method+'/'+basis+'/freq'+" GOODBYE!!"
        sys.exit(1)
    energies = extr_energy(out)  # calling the energy function
    freq_list = extr_entropy(out)
    S_sum = 0.0
    B_av = Bav(out)
    Sv_sum = 0.0
    Sr_sum = 0.0
    Sr_list = []
    Sv_list = []
    S_list = []
    for freq in freq_list:
        if float(freq) > 0.00:
            #print "\n\nThe frequency  "+str(freq)+": "
            I = mu(freq)
            mu_p = eff_mu(I,B_av)
            S_v=Sv(freq)
            Sv_list.append(S_v)
            Sv_sum = Sv_sum + S_v
            #print S_v
            S_r= Sr(mu_p)
            Sr_list.append(S_r)
            #print S_r
            Sr_sum = Sr_sum + S_r
            #print str(w(freq))+" * "+str(S_v)+" + " +str(1.0-w(freq))+" * "+str(S_r)
            S = w(freq)*S_v + (1.0-w(freq))*S_r
            S_list.append(S)
            S_sum = S_sum + S
        else:
            continue
    Final_entropy = S_sum*T
    Final_entropy_NOI = Sv_sum*T
    print '\n--------------------------------------'
    print "Interpolated: "+str(Final_entropy)+"kcal/mol"
    print '--------------------------------------'
    print "Non interpolated: "+str(Sv_sum*T) + " kcal/mol\n"
    #print '--------------------------------------'
    #print str(Sr_sum*T)

    #plt.plot(freq_list, S_list, 'g-')
    #plt.plot(freq_list, Sv_list, 'r-')
    #plt.plot(freq_list, Sr_list, 'b-')
    #plt.xlim([0, 350])
    #plt.show()

    #sys.exit(1)
    results1[d]=[] 
    results2[d]=[] 
    energies["S_vib_inter"] = Final_entropy/float(conv) 
    #energies["S_vib_NOI"] = Final_entropy_NOI/float(conv)
    energies["H_correct"] =  float(energies["Therm_correct"]) + float(energies["ThermH_correct"])
    energies['H_TS'] = energies["H_correct"] - energies["S_vib_inter"]
    if not d=='cat':
        results1[d]=[] 
        results2[d]=[] 
    for m in energies.iterkeys():
        if d == 'cat':
            cat_en[m] = energies[m]
            continue
        if d == "r":
           if ref:
               react_en[m] = energies[m]
               print "Extracted the reactant energy!"
               rel_en2 = 0.00
               results2[d].append([m, rel_en2])
               continue
           else: 
               pass
        if d == 'mon':
            mon_ener = energies[m]
            cat_mon_en[m] = float(mon_ener)+float(cat_en[m])
            energies[m] = cat_mon_en[m]
            #print m+" Energy of catalyst+monomer = "+str(cat_mon_en[m])
            if not ref:
                rel_en_cat_mon = 0.00
                results2[d].append([m,rel_en_cat_mon])
                continue
            else:
                pass

        if ref:
            rel_en2 = round((float(energies[m]) -float(react_en[m]))*float(conv) ,2)
        else:
            print "Entropy"+str(energies["S_vib_inter"])
            rel_en2 = round((float(energies[m])-float(cat_mon_en[m]))*float(conv) ,2)
        results2[d].append([m, rel_en2])
        print "******************"
        print '{:15s}  {:^55.2f}'.format(d+'('+m+')' ,rel_en2 )
    os.chdir(runDir)

print '\n\n'

#f = open('./data/energies_'+method+'_'+basis+'.dat2','w')
f = open('./data/thermo_data.dat','w')

f.write("#Method:  "+method+'\n')
f.write("#Basis:   "+basis+'\n' )
f.write('{:15s}'.format('Method   '))
for s in sorted(results1, key=molpro.struct_key):
    if 'mon' in s:
        pr = "Cat+mon"
        f.write('{:^10s}'.format(pr))
    elif 'cat' in s:
        continue
    else:
        f.write('{:^10.10s}'.format(s))
f.write('\n')

#for m in energies.iterkeys():
#    f.write('{:15s}'.format(m))
#    for s in sorted(results1, key=molpro.struct_key):
#        for i in results1[s]:
#            if i[0] == m:
#                energy = i[1]
#                if "." in str(energy):
#                    f.write('{:^10.2f}'.format(energy))
#                else:
#                    f.write('{:^10s}'.format(energy))
#    f.write('\n')
               
f.write('\n\n')
for m in energies.iterkeys():
    f.write('{:15s}'.format(m))
    for s in sorted(results2, key=molpro.struct_key):
        for i in results2[s]:
            if i[0] == m:
                energy = i[1]
                if "." in str(energy):
                    f.write('{:^10.2f}'.format(float(energy)))
                else:
                    f.write('{:^10s}'.format(float(energy)))
    f.write('\n')
               
f.write('\n\n')
f.close()


