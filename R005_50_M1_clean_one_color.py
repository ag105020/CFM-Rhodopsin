'''
Created on Apr 7, 2022

@author: keiin
'''


from pylab import *
from FigSetting2 import *
from sf import *
from st_data import *



rcParams.update({'lines.markersize':3,'lines.markeredgewidth':0})

def gft(a):
    return genfromtxt('../data/'+str(a)+'.csv',delimiter=',')

lat = gft('lat_fla')
surf = gft('M1_surf')
lay1 = gft('M1_lay1')
lay2 = gft('M1_lay2')
lay3 = gft('M1_lay3')
bot  = gft('M1_bot')

#OOOOOOOOOOOOOOOOOOOOOOOO
# Obtaining output
#OOOOOOOOOOOOOOOOOOOOOOOO

figure(1)
#plot(lat,bot,'o')
plot(lat,lay3,'o',color='#1F77B4')
plot(lat,lay2,'o',color='#1F77B4')
plot(lat,lay1,'o',color='#1F77B4')
plot(lat,surf,'o',color='#1F77B4')

ylabel('Rhodopsin per 100 genomes')
xlabel('latitude')
ylim(0,120)
title('Surface', y=1.02)

#========data part============
data = genfromtxt("../data/RhodopsinData_depth_separate.csv",delimiter=',').T

depth = data[0]
lat = data[1]
lon = data[2]
rhodopsin = data[3]/10

rcParams.update({'lines.markersize':10,
                 'lines.markeredgewidth':1})

plot(lat,rhodopsin,'o',color='k') 
#=============================

xticks((-50,-25,0,25,50,75))
xticklabels = ["50\N{DEGREE SIGN}S","25\N{DEGREE SIGN}S","0\N{DEGREE SIGN}","25\N{DEGREE SIGN}N","50\N{DEGREE SIGN}N","75\N{DEGREE SIGN}N"]
gca().set_xticklabels(xticklabels)

sf('Model1_One_color')

show()