'''
Created on Apr 7, 2022

@author: keiin
'''


from pylab import *
from FigSetting2 import *
from sf import *
from st_data import *
from E_BiosynthesisFromNH4 import *
from DG002_04_DOC import getDOC
from light import *


I = genfromtxt('../data/par_fla.csv',delimiter=',')
kd = genfromtxt('../data/kd_fla.csv',delimiter=',')
lat = genfromtxt('../data/lat_fla.csv',delimiter=',')
lon = genfromtxt('../data/lon_fla.csv',delimiter=',')
DOC = genfromtxt('../data/DOC_fla.csv',delimiter=',')

I = I*1e6/86400 #converting from mol m-2 d-1 to umol m-2 s
z = 0
I = I*exp(-kd*z)
#Light coefficient
AI = 0.00863364097132997  #(umol m-2 s-1) Light coefficient

m = 0.08
m2 = 0.1

#Tunable parameters
Rmax = 9.226909902953164 #(d-1) Maximum energy gain by rhodopsin
aR = 0.006961874073858994 #(d-1/uM) Affinity
a =  0.01801081555475789 #(d-1/uM) Affinity
ep = 0.2881772139121359  #McCarty 2007
E = E_BiosynthesisFromNH4(ep) #(mol C mol C-1) CO2 production rate

#Rhodopsin bacteria
R = Rmax*(1-exp(-AI*I)) #(d-1) energy gain by rhodopsin
R[isnan(DOC)==1] = nan
VcR = aR*DOC  #(d-1) DOC uptake rate
print(min(VcR))
MuR = (VcR+R)/(1+E)

c = R>MuR*E #c: condition
MuR[c] = VcR[c]

#Non-rhodopsin bacteria
Vc = a*DOC #(d-1) DOC uptake rate
Mu = Vc/(1+E) #(d-1) Growth rate
Ec = Mu*E #(d-1) Eneryg consumption

#OOOOOOOOOOOOOOOOOOOOOOO
#Ecosystem model
#OOOOOOOOOOOOOOOOOOOOOOO
XR = (MuR - m)/m2
X = (Mu - m)/m2
X[X<0] = nan
XR[XR<0] = nan
fXR = XR/(XR + X) #(dimensionless) fraction of XR
fXR[fXR>1] = nan
fXR[fXR<0] = nan
fXRpercent = fXR*100



#OOOOOOOOOOOOOOOOOOOOOOOO
# Obtaining output
#OOOOOOOOOOOOOOOOOOOOOOOO

figure(1)
plot(lat,fXRpercent,'o',)
ylabel('Rhodopsin %')
xlabel('latitude')
ylim(0,100)
title('Model1')

sf('Model1C')
st(fXRpercent,'M1_surf')

show()