#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:02:24 2022
@author: stillwel
"""

import numpy

def PoissThin(a):
    s = numpy.array(a)
    s.astype(numpy.int64,casting='unsafe')
    f = numpy.random.binomial(s, 0.5, size = s.shape)
    return(f)

if __name__ == '__main__':
    try:
        x
    except:
        x = [10,10,10]
    F = PoissThin(x)