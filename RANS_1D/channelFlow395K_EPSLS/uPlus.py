# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 17:40:55 2018

@author: constantinos
"""

"All combined Plots of the high Reynolds number models "

##producing plots to compare DNS DATA
#must be in same folder as the files
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os 
#plt.rcParams.update({'font.size': 14})
#allows to pass arguments when calling the script
#import argparse
import sys

#i.e.  python3 scriptName.py firstArgument secondArgument....  nArguments
#sys.argv[1] = firstArgument
#sys.argv[2] = secondArgument
#sys.argv[n] = secondArgument
#def main(argv):

#try:
#  runDirectory = sys.argv[1]
#except Exception:
#  print("You did not specify a file but I will plot the DNS")
#  pass
  #sys.exit(1)
  #Vertical channel

def main(argv):
  runDirectory = argv
  def utau(dUdy,nu):
      tw = nu*dUdy #from low-Re
      u_tau = math.sqrt(tw/1)
      print(u_tau)
      return u_tau
      

  def yplus(ut,nu,y ):
     yplus = (y*ut)/nu
     return yplus

  def Uplus(ut,u):
     return u/ut
          
  variableColumnIndex =[8]#,10,5,2]  
  dnsVariableColumnIndex =[2]#,9,6,8]        
  yAxislabelsList =['$U^{+}$']#,'$T^{+}$','-$\overline{uv}^{+}$','$k^{+}$' ]           
  modelList = ['$k-{\epsilon}$ SGDH']           
  markerColorList=['r-x','y+-', 'gv-','m*-']
      #markerColorList=['r--','y-', 'gv-','m*-']
  dnsdata =pd.read_excel('/home/mariawsl/CFD_scratch/OpenFoam/channelFlow/DNS_Retau_395/Poi395_4th_Cuhf071.xls',comment='#')
          
  #Plot figure
  fig = plt.figure(figsize=(9.5, 9.5))
  fig.suptitle('$Re_{tau}$=395, $Re_b$~14147',  fontsize =16)
  fig.subplots_adjust(top=0.9)
  fig.subplots_adjust(wspace=0.3, hspace = 0.25)
       
  for varIndex in range(len(variableColumnIndex)):
           
     #axes1 = fig.add_subplot(2,2,varIndex+1)
     axes1 = fig.add_subplot(2,2,1)
     #DNS data
     DNSyplus = dnsdata.iloc[1:224,1].values
     DNSvariable = dnsdata.iloc[1:224,2].values
     #print(DNSvariable)
     #DNSvariable = dnsdata.iloc[1:,dnsVariableColumnIndex[varIndex]].values
     #DNSuv = dnsdata.iloc[1:224,4].values
     #DNSTKE = dnsdata.iloc[1:224,8].values
     #DNSTplus = dnsdata.iloc[1:224,3].values
     axes1.plot(DNSyplus,  DNSvariable, 'ko', label='DNS Kawamura',linewidth=4)
     #Numerical simulation data    
     for modelIndex in range(len(modelList)):  
     #dataRead = pd.read_table('profile_code.dat',sep='\s+')
         #dataRead = pd.read_csv(runDirectory,skiprows=1)
         dataRead = pd.read_csv(runDirectory)
         U1       = dataRead.iloc[0,1]#.values #should be zero
         U2	      = dataRead.iloc[1,1]
         y1       = dataRead.iloc[0,0]
         y2       = dataRead.iloc[1,0]
         
         dUdy     = (U2-U1)/(y2-y1)
         print(dUdy)
         nu       = 0.002531
         length   =int(dataRead[dataRead.columns[0]].count())#count rows and convert to int
                  
     # if varIndex==1:
     #    variableModel = 34.7 -(dataRead.iloc[0:16,variableColumnIndex[varIndex]].values)
         Uplus         = Uplus(utau(dUdy,nu),dataRead.iloc[1:90,1].values)     
         yplusModel    = yplus(utau(dUdy,nu),nu, dataRead.iloc[1:90,0].values)
         axes1.plot(yplusModel,Uplus,markerColorList[modelIndex],label=modelList[modelIndex], linewidth=3,ms=5)
            
     #if varIndex < 2:
     axes1.set_xscale('log')#change the x-scale to logarithmic
     axes1.set_xlabel('$y^{+}$',fontsize =16)
     axes1.set_ylabel(yAxislabelsList[varIndex],fontsize =16)        
             
     handles,labels =axes1.get_legend_handles_labels()
         #axes4.legend(bbox_to_anchor=(1.22,2), shadow=True, fontsize = 10)
     lgd=axes1.legend(handles,labels,loc='upper center',bbox_to_anchor=(0.5,-0.2), shadow=True, fontsize = 12)

     #plt.show()    
     plt.savefig('Retau395Uplus.png',dpi=150,bbox_extra_artists=(lgd,),bbox_inches='tight')
#try:
# if __name__=="__main__":
#  main(sys.argv[1])
#  print("Arg was given, I will proceed")
#except Exception:
#  print("You did not specify a file and I will not plot anything")
##  pass
#  sys.exit(1)
    #Vertical channel
if __name__ == "__main__":
    main(sys.argv[1])
