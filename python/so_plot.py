#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 11:59:18 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

fd = open('percent_left_code13.dat', 'rb')
fd = open('so_tot_code13.dat', 'rb')
#fd = open('o23_tot_code13.dat', 'rb')
#fd = open('o21_tot_code13.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double)
fd.close()
so=datain #np.flip(datain)



time=np.arange(600)

  

plt.plot(time,so)
plt.show()