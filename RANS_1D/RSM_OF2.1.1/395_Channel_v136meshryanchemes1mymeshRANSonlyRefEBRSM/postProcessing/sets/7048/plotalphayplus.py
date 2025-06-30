#!/usr/bin/python3.5
import matplotlib.pyplot as plt
import numpy as np

#prepare file
target = open('leftPatch_alpha_prep.xy', 'w')
with open('leftPatch_alpha.xy', 'r') as f:
    for line in f:
         line = line.split()
         target.write(line[0])
         target.write('   ')
         target.write(line[1])
         target.write('\n'),
target.close()


#plot file
y,alpha=(np.loadtxt('leftPatch_alpha_prep.xy',delimiter='   ',unpack=True))
#w,z=(np.loadtxt('Poi395_4th_D_uhf071uithetaprep.dat',delimiter='   ',unpack=True))
#w=w*392.24
yplus=y*(5.84404e-05**0.5)/0.00002
point96loc=0.180945
point96locplus=point96loc*(5.84404e-05**0.5)/0.00002

plt.figure(1)
plt.plot(yplus,alpha,label='RANS')
#plt.plot(point96locplus,alpha,label='0.96')
#plt.axvline(x=point96locplus)
#setlimits,locationimportant
#plt.xaxis([0, 1])
#plt.xlim(0, 1)
#plt.ylim(0, 0.1) 
plt.yticks(np.arange(min(y), max(y)+0.1, 0.1))
plt.grid()
plt.xlabel('y (m)')
plt.ylabel('T+')
plt.title('T+ vs y\n')
plt.legend()
#locationofsavefigbeforeshowisimportant
plt.savefig("vtheta.png")
plt.show()

plt.figure(2)
plt.semilogx(yplus,alpha,label='RANS')
plt.yticks(np.arange(min(y), max(y)+0.1, 0.1))
#plt.axvline(x=point96locplus)
plt.grid()
plt.xlabel('y+')
plt.ylabel('alpha')
plt.title('alpha vs y+\n')
plt.savefig("vthetaflog.png")
plt.show()
 