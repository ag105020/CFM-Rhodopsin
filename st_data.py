'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def st(parameter,name):
    First_part="..\\data\\"
    Second_part=str(name)
    Last_part=".csv"
    savetxt(First_part+Second_part+Last_part, parameter, delimiter=",",fmt='%.8e')