'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def st(parameter,name):
    First_part="C:\\Users\\Keiin\\Desktop\\figures\\02\\19 Rhodopsin bacteria\\"
    Second_part=str(name)
    Last_part=".csv"
    savetxt(First_part+Second_part+Last_part, parameter, delimiter=",",fmt='%.8e')