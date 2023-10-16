'''
Created on Jan 20, 2022

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from E_BiosynthesisFromNH4 import *
from sf import *
from st import *

#Data=======
data = genfromtxt("../data/RhodopsinData.csv",delimiter=',').T
depth = data[0]
lat = data[1]
lon = data[2]
rhodopsin = data[3]

#Data=======
z = genfromtxt("../data/z.csv",delimiter=',').T[1,:13]   #(m) Depth
DOC = genfromtxt("../data/DOCz.csv",delimiter=',')[:13]  #(uM) DOC concentration for z direction

#Light profile: From Masuda et al 2020
z0 = 30         #(m) Depth of 1/e light of surface
I0 = 340.51       #(umol m-2 s-1) Surface light intensity
I = I0*exp(-z/z0) #(umol m-2 s-1)
AI = 0.00863364097132997  #(umol m-2 s-1) Light coefficient

m = 0.08
m2 = 0.1

#Tunable parameters
Rmax = 9.226909902953164 #(d-1) Maximum energy gain by rhodopsin
aR = 0.006961874073858994 #(d-1/uM) Affinity
a = 0.01801081555475789 #(d-1/uM) Affinity
ep = 0.2881772139121359  #McCarty 2007
E = E_BiosynthesisFromNH4(ep) #(mol C mol C-1) CO2 production rate

aR1max = 0.008123184721478018  #(d-1/uM) Increased affinity with light (for rhodopsin prokaryote)

aCfixMax = 0.31/0.69 #(dimensionless) 31% of C fixation (Palovaara et al 2014)

def function(aR1max,aCfixMax,Rmax):
    aR1 = aR1max*(1-exp(-AI*I))  
    aCfix = aCfixMax*(1-exp(-AI*I))
    #Rhodopsin bacteria
    R = Rmax*(1-exp(-AI*I)) #(d-1) energy gain by rhodopsin
    VcR = (aR+aR1)*DOC*(1+aCfix)  #(d-1) DOC uptake rate
    
    MuR = (VcR+R)/(1+E)
    
    c = R>MuR*E #c: condition
    MuR[c] = VcR[c]
    
    EcR = MuR*E #(d-1) Energy consumption for Rhodopsin bacteria
    R[c] = MuR[c]*E
    EcR[c] = EcR[c] - R[c]
    
    which = zeros(size(z))
    which[c] = 1
    
    #Non-rhodopsin bacteria
    Vc = a*DOC #(d-1) DOC uptake rate
    Mu = Vc/(1+E) #(d-1) Growth rate
    Ec = Mu*E #(d-1) Energy consumption
    
    #OOOOOOOOOOOOOOOOOOOOOOO
    #Ecosystem model
    #OOOOOOOOOOOOOOOOOOOOOOO
    XR = (MuR - m)/m2
    X = (Mu - m)/m2
    fXR = XR/(XR + X) #(dimensionless) fraction of XR
    
    return fXR, MuR, Mu, Ec, EcR, R



fXR3, MuR3, Mu3, Ec3, EcR3, R3 = function(aR1max,aCfixMax,Rmax)
fXR4, MuR4, Mu4, Ec4, EcR4, R4 = function(aR1max,aCfixMax,Rmax*0.8)

#OOOOOOOOOOOOOOOOOOOO
# Plotting
#OOOOOOOOOOOOOOOOOOOO

mean_rhodopsin = genfromtxt("../data/Mean_Rhodopsin.csv",delimiter=',')

figure(2)

plot(fXR3*100,z,"-",label='Default $\mathit{R_{max}}$',color="#FF7F0E")
plot(fXR4*100,z,"-.",label='20% lower $\mathit{R_{max}}$',color="k")
ylabel('Depth (m)')
xlabel('Rhodopsin per 100 genomes')
legend(edgecolor='k',loc=4)
xlim(-3,120)
gca().invert_yaxis()
sf('Sensitivity')


show()
