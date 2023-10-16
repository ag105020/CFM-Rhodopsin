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
Bottop = gft('M2_Bot_top')
lay6 = gft('M2_lay6')
lay7 = gft('M2_lay7')
lay8 = gft('M2_lay8')
lay9 = gft('M2_lay9')
lay10 = gft('M2_lay10')
lay11 = gft('M2_lay11')
lay12 = gft('M2_lay12')



#OOOOOOOOOOOOOOOOOOOOOOOO
# Obtaining output
#OOOOOOOOOOOOOOOOOOOOOOOO

figure(1)
#plot(lat,bot,'o')
plot(lat,lay12,'o',color='#1F77B4')
plot(lat,lay11,'o',color='#1F77B4')
plot(lat,lay10,'o',color='#1F77B4')
plot(lat,lay9,'o',color='#1F77B4')
plot(lat,lay8,'o',color='#1F77B4')
plot(lat,lay7,'o',color='#1F77B4')
plot(lat,lay6,'o',color='#1F77B4',label='Model')
#plot(lat,Bottop,'o',color='#1F77B4')

ylabel('Rhodopsin per 100 genomes')
xlabel('latitude')
ylim(0,120)
title('Deep', y=1.02)

#========data part============
data = genfromtxt("../data/RhodopsinData_depth_separate.csv",delimiter=',').T

depth = data[0]
lat = data[9]
lon = data[10]
rhodopsin = data[11]/10

rcParams.update({'lines.markersize':10,
                 'lines.markeredgewidth':1})

plot(lat,rhodopsin,'o',color='k',label='Data') 
#=============================
legend(edgecolor='k')
xticks((-50,-25,0,25,50,75))
xticklabels = ["50\N{DEGREE SIGN}S","25\N{DEGREE SIGN}S","0\N{DEGREE SIGN}","25\N{DEGREE SIGN}N","50\N{DEGREE SIGN}N","75\N{DEGREE SIGN}N"]
gca().set_xticklabels(xticklabels)

sf('Model2_bot_lat2')

show()