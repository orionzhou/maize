#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re

from maize.apps.base import eprint
from itertools import chain
flatten = chain.from_iterable

LEFT, RIGHT = 1, -1
ptn_loc = re.compile("^([\w\-\.]+)\:([\d,]+)(\-|\.{1,2})([\d,]+)$")

def locStr2Ary(locS):
    locS = str(locS)
    if not locS:
        return []
    locA = []
    try:
        for loc in locS.split(","):
            beg, end = loc.split("-")
            locA.append([int(beg), int(end)])
        return locA
    except:
        print("invalid locStr: %s" % locS)

def locAry2Str(locA):
    if not locA or len(locA) == 0:
        return ""
    try:
        return ",".join(["%s-%s" % (x[0], x[1]) for x in locA])
    except:
        print("invalid locAry: ", locA)

def locAryLen(locA):
    if not locA or len(locA) == 0:
        return 0
    try:
        return sum([x[1] - x[0] for x in locA])
        #return sum([x[1] - x[0] + 1 for x in locA])
    except:
         print("invalid locAry: ", locA)

def join_ranges(data, offset=0):
    data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT)) \
            for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0: 
            yield x, value - offset


if __name__ == '__main__':
    print(locStr2Ary("7-10,8-19"))
    locStr2Ary(7)
    print(locAry2Str([[8,9], [10,11]]))
    locAry2Str([[8,9], [11]])
    print(locAryLen([[8,9], [10,11]]))
