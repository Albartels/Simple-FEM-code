#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# This function gives an interpolated value back for the discretized loadcurve
#
# coded by abartels
###############################################################################

import numpy as np

def scale_factor(loadcurve, time):
    
    startvalue=0.0
    slope=0.0
    
    for ii in range(np.shape(loadcurve)[1]-1):
        if (time>=loadcurve[0][ii]) and (time<loadcurve[0][ii+1]):
            slope=(loadcurve[1][ii+1]-loadcurve[1][ii])/\
                (loadcurve[0][ii+1]-loadcurve[0][ii])
            startvalue=loadcurve[1][ii]
    return (startvalue+slope*time)