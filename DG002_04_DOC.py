'''
Created on Apr 6, 2022

@author: keiin
'''


from pylab import *
from scipy.io.netcdf import netcdf_file

def getDOC():
    nc = netcdf_file('../data/DOC_Lonborg_etal2018_K.nc')
    DOC = nc.variables['DOC_R1'][:12,:,:]
    
    return DOC

# 
# DOC = getDOC()
# savetxt('DOC.csv',DOC[0],delimiter=',',fmt='%.8e')
