# -*- coding: utf-8 -*-
""""
Created on Mon Aug 29 11:42:00 2022

@author: constantinos
"""

###################import modules#############################################
from pyexpat.errors import XML_ERROR_DUPLICATE_ATTRIBUTE
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import math
import numpy as np
#############################################################################

#class to create each graph
class createGraph:
    #constructor and attributes of objects
    def __init__(self,title,xAxisLbl):
        #return
        self.title  = title
        self.xAxisLbl=xAxisLbl
            
    #member function to plot graph
    def plotGraph(self,xdata,yplus):

        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(top=0.96, left=0.2,bottom=0.15,right=0.7) 

        plt.plot(xdata,yplus)
        plt.xlabel(self.xAxisLbl,fontsize =16)
        plt.savefig(str(self.title)+'.png',dpi=100,bbox_inches="tight")
        

pathU='/home/mariawsl/CFD_scratch/OpenFoam/channelFlow/channelFlow395K_EPSLS/postProcessing/sample/100'
file = '/verticalline_U.csv'

readFile=pd.read_csv(pathU+file,skiprows=1)
U       = readFile.iloc[:,1]
y       = readFile.iloc[:,0]
varTitle= "Retau_395"
labelX  = "U"


plot1=createGraph(varTitle,labelX)
plot1.plotGraph(U,y)

#if __main__ == "__main__":
#    main(sys.argv[1])