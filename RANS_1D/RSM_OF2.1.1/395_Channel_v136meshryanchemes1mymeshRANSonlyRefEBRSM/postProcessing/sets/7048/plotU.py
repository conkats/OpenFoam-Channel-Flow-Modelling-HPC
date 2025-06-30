#!/usr/bin/python3.5
import matplotlib.pyplot as plt
import numpy as np

#prepare file
target = open('leftPatch_Uprep.xy', 'w')
with open('leftPatch_U.xy', 'r') as f:
    for line in f:
         line = line.split()
         target.write(line[0])
         target.write('   ')
         target.write(line[1])
         target.write('\n'),
target.close()

target = open('Poi395_4th_DUprep.dat', 'w')
with open('Poi395_4th_D.dat', 'r') as f:
     for line in f:
          line = line.split()
          target.write(line[1])
          target.write('   ')
          target.write(line[2])
          target.write('\n'),
target.close()

#plot file
WSS=0.0000584404
x,y=(np.loadtxt('leftPatch_Uprep.xy',delimiter='   ',unpack=True))
w,z=(np.loadtxt('Poi395_4th_DUprep.dat',delimiter='   ',unpack=True))
#w=w*392.24
y=y/(WSS**(0.5))
n=WSS**(0.5)/(0.00002)
x=x*n


plt.figure(1)
plt.plot(x,y,label='Hybrid RANS',color='b')

plt.plot(w,z,label='DNS',linestyle=' ',marker='*',color='k')
#setlimits,locationimportant
#plt.xaxis([0, 1])
#plt.xlim(0, 1)
#plt.ylim(0, 0.1) 
plt.grid()
plt.xlabel('y+')
plt.ylabel('U+')
plt.title('EBRSM U+ vs y+\n')
#plt.legend()
#locationofsavefigbeforeshowisimportant
plt.savefig("Uf.png")
plt.show()

plt.figure(2)
plt.semilogx(x,y,label='Hybrid RANS',color='b')
plt.semilogx(w,z,label='DNS',linestyle=' ',marker='*',color='k')
plt.grid()
plt.xlabel('log y+')
plt.ylabel('U+')
plt.title('EBRSM U+ vs log y+\n')
plt.savefig("Uflog.png")
plt.show()
 