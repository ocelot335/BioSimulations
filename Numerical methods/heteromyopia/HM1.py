import os
import numpy
import string
import pandas as pd

data = pd.read_csv('dset1_all_in_one_HM.dat.csv')
N1_sim = data['N 1']
N2_sim = data['N 2']

min_error = float('inf')
best_al = None
plot1 = open('plot_HM1.txt', 'a')
plot1.truncate(0)

for al in numpy.linspace(0.4, 0.4, num=1):
    error = 0
    correct_ratios = []
    i = 0
    for sd_in in numpy.linspace(0.01, 0.22, num=8):
        for sd_ex in numpy.linspace(0.01, 0.22, num=8):
            cmd = f'./a point -t test1 -dim 1 -al {al} -points 1024 -a 10.000000 -b 0.4 0.4 -dvec 0.2 0.2 -dmat 0.001 0.001 0.001 0.001 -sw {sd_in} {sd_ex} {sd_ex} {sd_in} -sm 0.04 0.06'
            os.system(cmd)
            n = open('test1.N', 'r')
            s = (n.readline()).split(' ')
            print("-----")
            print(str(sd_in)+' '+str(sd_ex))
            N1 = 0 if s[0] == '-nan' else float(s[0])
            N2 = 0 if s[1] == '-nan' else float(s[1])
            N1 = 0 if N1 < 0 else N1
            N2 = 0 if N2 < 0 else N2
            plot1.write(str(sd_in) + ' ' + str(sd_ex) + ' ')
            plot1.write(str(10*N1)+' ')
            plot1.write(str(10*N2)+'\n')
            n.close()

            print(str(10*N1)+' '+str(10*N2)+' '+str(N1_sim[i])+' '+str(N2_sim[i])+' ')
            error += (N1*10 - N1_sim[i])**2 + (N2*10 - N2_sim[i])**2

            if abs(N1*10 - N1_sim[i]) <= max(100, 0.1*N1_sim[i]) and abs(N2*10 - N2_sim[i]) <= max(100, 0.1*N2_sim[i]):
                correct_ratios.append(1)
            else:
                correct_ratios.append(0)

            i += 1

    if error < min_error:
        min_error = error
        best_al = al

    with open('errors_al_HM1.txt', 'a') as f:
        f.write(str(al)+' ')
        f.write(str(error))
        f.write('\n')

    with open(f'correct_ratios_al_HM1.txt', 'a') as f:
        f.write(str(al)+' ')
        f.write(str(sum(correct_ratios) / len(correct_ratios)))
        f.write('\n')

print(f'Лучшее значение параметра -al: {best_al}')
