'''
Created on Jan 20, 2022

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from E_BiosynthesisFromNH4 import *
from sf import *
from st import *


#OOOOOOOOOOOOOOOOOOOO
# Plotting
#OOOOOOOOOOOOOOOOOOOO

PImax = 68
a = 0.5
I = arange(0,500,0.1)
PI = PImax*(1-exp(-a*I))


mean_rhodopsin = genfromtxt("../data/Mean_Rhodopsin.csv",delimiter=',')


I_data = genfromtxt("../data/Rhodopsin_Light.csv",delimiter=',').T
print(I_data)
figure(9)
plot(I_data[2],I_data[3],'o',markersize=8,label='Data')
plot(I,PI,label="PI")
xlabel('Light intensity ($\mu$mol m$^{-2}$ s$^{-1}$)')
ylabel('Rhodopsin per 100 genomes')
legend(edgecolor='k')
sf('PI_fit')


show()
