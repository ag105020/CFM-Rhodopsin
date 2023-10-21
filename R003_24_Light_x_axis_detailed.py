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
z = genfromtxt("../data/z2.csv",delimiter=',').T[1,:13]   #(m) Depth
DOC = genfromtxt("../data/DOCz2.csv",delimiter=',')[:13]  #(uM) DOC concentration for z direction

#Light profile: From Masuda et al 2020
z0 = 30         #(m) Depth of 1/e light of surface
I0 = 340.51       #(umol m-2 s-1) Surface light intensity
I = I0*exp(-z/z0) #(umol m-2 s-1)
AI = 0.00863364097132997  #(umol m-2 s-1) Light coefficient
I[0] = 480
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

def function(aR1max,aCfixMax):
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


fXR1, MuR1, Mu1, Ec1, EcR1, R1 = function(0,0)
fXR2, MuR2, Mu2, Ec2, EcR2, R2 = function(aR1max,0)
fXR3, MuR3, Mu3, Ec3, EcR3, R3 = function(aR1max,aCfixMax)

#OOOOOOOOOOOOOOOOOOOO
# Plotting
#OOOOOOOOOOOOOOOOOOOO

# figure(1)
# plot(XR,z,label='Rhodopsin')
# plot(X,z,label='Non-Rhodopsin')
# xlabel('Population')
# xscale('log')
# ylabel('Depth (m)')
# legend(edgecolor='k')
# gca().invert_yaxis()
# sf('Population')

mean_rhodopsin = genfromtxt("../data/Mean_Rhodopsin.csv",delimiter=',')

figure(2)
plot(rhodopsin/10,depth,'o',markersize=8,label='Data')
plot(mean_rhodopsin/10,z[:-1],'ro',label='Mean Data')
plot(fXR1*100,z,"--",label='Model1',color="#FF7F0E")
plot(fXR2*100,z,label='Model2',color="#FF7F0E")
plot(fXR3*100,z,":",label='Model3',color="#FF7F0E")
ylabel('Depth (m)')
xlabel('Rhodopsin per 100 genomes')
legend(edgecolor='k')
gca().invert_yaxis()
#sf('Model_data')

# print(fXR*100)
# print(z[:-1])
# st(z[:-1],'depth')
# st(fXR*100,'Rhodo_ratio')
#
# figure(3)
# plot(MuR,z,label='Rhodopsin')
# plot(Mu,z,label='Non-Rhodopsin')
# ylabel("Depth (m)")
# xlabel('$\mu$ (d$^{-1}$)')
# legend(edgecolor='k')
# gca().invert_yaxis()
# sf('Growth rate')

figure(41,figsize=(6.5,8))
stackplot(z,MuR1,EcR1,R1)
xlabel("Depth (m)",rotation=180)
ylabel("C allocation (d$^{-1}$)")
ylim(top=4)
xticks(rotation=90)
yticks(rotation=90)
#sf('C allocation Rhodopsin M1')

figure(42,figsize=(6.5,8))
stackplot(z,MuR2,EcR2,R2)
xlabel("Depth (m)",rotation=180)
ylabel("C allocation (d$^{-1}$)")
ylim(top=4)
xticks(rotation=90)
yticks(rotation=90)
#sf('C allocation Rhodopsin M2')

figure(43,figsize=(6.5,8))
stackplot(z,MuR3,EcR3,R3)
xlabel("Depth (m)",rotation=180)
ylabel("C allocation (d$^{-1}$)")
ylim(top=4)
xticks(rotation=90)
yticks(rotation=90)
#sf('C allocation Rhodopsin M3')
#
#
figure(51,figsize=(6.5,8))
stackplot(z,Mu1,Ec1)
xlabel("Depth (m)",rotation=180)
ylabel("C allocation (d$^{-1}$)")
ylim(top=4)
xticks(rotation=90)
yticks(rotation=90)
#sf('C allocation Non-Rhodopsin M1')

# figure(52,figsize=(6.5,8))
# stackplot(z,Mu2,Ec2)
# xlabel("Depth (m)",rotation=180)
# ylabel("C allocation (d$^{-1}$)")
# ylim(top=4)
# xticks(rotation=90)
# yticks(rotation=90)
#
# sf('C allocation Non-Rhodopsin M2')
#
# figure(53,figsize=(6.5,8))
# stackplot(z,Mu3,Ec3)
# xlabel("Depth (m)",rotation=180)
# ylabel("C allocation (d$^{-1}$)")
# ylim(top=4)
# xticks(rotation=90)
# yticks(rotation=90)
#sf('C allocation Non-Rhodopsin M3')

# figure(6)
# Colors = ['#1F77B4','#FF7F0E','#2CA02C']
# Names = ['Biomass','Respiration','C saving by Rho.']
# stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
# legend(loc=1,edgecolor = 'k')
# sf('legend')
#
# figure(7)
# plot(VcR[c],z[c])
# gca().invert_yaxis()

figure(8)
plot(MuR1/Mu1,z,"--",color="#1F77B4",label='Model1')
plot(MuR2/Mu2,z,color="#1F77B4",label='Model2')
plot(MuR3/Mu3,z,":",color="#1F77B4",label='Model3')
ylabel("Depth (m)")
xlabel('$\mu_{R}$/$\mu_{NR}$ ')
legend(edgecolor='k')
gca().invert_yaxis()
sf('Growth rate ratio')


I_data = genfromtxt("../data/Rhodopsin_light_C.csv",delimiter=',').T
print(I_data)
figure(9)
plot(I_data[0]*1e6/86400,I_data[1],'o',markersize=8,label='Data')
plot(I, fXR1*100,"--",label='Model1',color="#FF7F0E")
plot(I, fXR2*100,label='Model2',color="#FF7F0E")
plot(I, fXR3*100,":",label='Model3',color="#FF7F0E")
xlabel('Light intensity ($\mu$mol m$^{-2}$ s$^{-1}$)')
ylabel('Rhodopsin per 100 genomes')
legend(edgecolor='k')
sf('Model_data (light)')

figure(10)
plot(I, MuR1/Mu1,"--",color="#1F77B4",label='Model1')
plot(I, MuR2/Mu2,color="#1F77B4",label='Model2')
plot(I, MuR3/Mu3,":",color="#1F77B4",label='Model3')
xlabel('Light intensity ($\mu$mol m$^{-2}$ s$^{-1}$)')
ylabel('$\mu_{R}$/$\mu_{NR}$ ')
legend(edgecolor='k')
ylim(top=3.8)
#gca().invert_yaxis()
sf('Growth rate ratio (light)')




show()
