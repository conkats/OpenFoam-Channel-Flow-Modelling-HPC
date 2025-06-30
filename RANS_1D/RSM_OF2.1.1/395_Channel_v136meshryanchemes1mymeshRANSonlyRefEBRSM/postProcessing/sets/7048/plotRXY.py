#!/usr/bin/python3.5
import matplotlib.pyplot as plt
import numpy as np

#prepare file
target = open('leftPatch_TotRXYprep.xy', 'w')
with open('leftPatch_R.xy', 'r') as f:
    for line in f:
         line = line.split()
         target.write(line[0])
         target.write('   ')
         target.write(line[2])
         target.write('\n'),
target.close()

target = open('uvprep.xy', 'w')
with open('uv', 'r') as f:
    for line in f:
         line = line.split()
         target.write(line[0])
         target.write('   ')
         target.write(line[1])
         target.write('\n'),
target.close()

target = open('Poi395_4th_Duvprep.dat', 'w')
with open('Poi395_4th_Dv2uv.dat', 'r') as f:
     for line in f:
          line = line.split()
          target.write(line[1])
          target.write('   ')
          target.write(line[3])
          target.write('\n'),
target.close()

#plot file
x3,y3=(np.loadtxt('leftPatch_TotRXYprep.xy',delimiter='   ',unpack=True))
x4,y4=(np.loadtxt('uvprep.xy',delimiter='   ',unpack=True))

w3,z3=(np.loadtxt('Poi395_4th_Duvprep.dat',delimiter='   ',unpack=True))

#x1,y1=(np.loadtxt('../../../Varattempt1/graphs/14000/TotVarprep.xy',delimiter='   ',unpack=True))
##w,z=(np.loadtxt('Poi395_4th_D_uhf071prep.dat',delimiter='   ',unpack=True))
#w=w*392.24
##z=z**(2)
WSS=0.0000584404
y3=y3/(WSS)*-1
n3=WSS**(0.5)/(0.00002)
x3=x3*n3

y4=y4/(0.000058822)*-1
n4=0.000058822**(0.5)/(0.00002)
x4=x4*n4
#y1=y1/(0.36298012967**(2))
#n1=0.0000567423**(0.5)/(0.00002)
#x1=x1*n1
plt.figure(1)
plt.plot(x3,y3,label='<uv>+',color='b')
plt.plot(x4,y4,label='<uv>+ StarCCM+',linestyle='--',color='b')
#plt.plot(x3,k,label='k')

plt.plot(w3,z3,label='<uv>+ DNS',linestyle=' ',marker='*', color='k')
##plt.plot(w,z,label='DNS',linestyle=' ',marker='*')
#plt.plot(x1,y1,label='TMean0')
#setlimits,locationimportant
#plt.xaxis([0, 1])
#plt.xlim(0, 1)
#plt.ylim(0, 0.1) 
plt.grid()
plt.xlabel('y+')
plt.ylabel('<uv>+')
plt.title('EBRSM <uv>+ vs y+\n')
plt.legend()
#locationofsavefigbeforeshowisimportant
plt.savefig("uvf.png")
plt.show()

plt.figure(2)
plt.semilogx(w,z,label='Rxx',linestyle=' ',marker='*')
plt.semilogx(w1,z1,label='Rzz',linestyle=' ',marker='*')
plt.semilogx(w2,z2,label='Ryy',linestyle=' ',marker='*')
##plt.semilogx(w3,z3,label='Rxy',linestyle=' ',marker='*')
##plt.semilogx(w,z,label='DNS',linestyle=' ',marker='*')
#plt.semilogx(x1,y1,label='TMean0')
plt.grid()
plt.xlabel('y+')
plt.ylabel('U+')
plt.title('Hybrid LES T+ vs y+\n')
plt.savefig("uvflog.png")
plt.show()
