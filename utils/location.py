#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re

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
        return sum([x[1] - x[0] + 1 for x in locA])
    except:
         print("invalid locAry: ", locA)

def maketile(beg, end, winsize, winstep):
    import math
    size = end - beg + 1
    nwin = math.ceil((size - winsize) / winstep) + 1
    
    mergelast = False
    if nwin > 1:
        size_last = end - (beg + winstep * (nwin-2) + winsize) + 1
        #print(size_last)
        if float(size_last) / winsize < 0.5:
            mergelast = True

    wins = []
    for i in range(0, nwin):
        wbeg = beg + winstep * i
        wend = beg + winstep * i + winsize - 1
        wend = min(end, wend)
        if i == nwin - 2 and mergelast:
            wend = end
        elif i == nwin - 1 and mergelast:
            continue
        wins.append([wbeg, wend])
    return wins

if __name__ == '__main__':
    print(locStr2Ary("7-10,8-19"))
    locStr2Ary(7)
    print(locAry2Str([[8,9], [10,11]]))
    locAry2Str([[8,9], [11]])
    print(locAryLen([[8,9], [10,11]]))
