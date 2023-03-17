#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:17:00 2021

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np

fd = open('flu_image.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(200,200)
fd.close()
data=np.flipud(datain)


plt.imshow(datain,cmap='Reds_r')


plt.show()
