'''
Created on Jul 13, 2022
Written by Meng Gao, revised by K
@author: 31417
'''
from pylab import *


def Ical():
    I0 = genfromtxt("../Data/LightData.csv",delimiter=',').T
    zn=12
    xn=360
    yn=180
    I = zeros([zn,xn,yn])
    z = genfromtxt("../Data/Depth.csv",delimiter=',')
    z0 = 30
    for i in range(size(z)):
        I[i] = I0 * exp(-z[i]/z0)
    return(I*1000000/86400)     

I = Ical()

# for i in range(12):
#     figure(i)
#     pcolormesh(arange(0,360),arange(-90,90),I[i].T,shading = 'flat')
#     colorbar()
# show()