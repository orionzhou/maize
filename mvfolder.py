#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'recursively move files/folders'
    )
    parser.add_argument(
            'dirw', nargs = '?', default = '/scratch.global/zhoux379', help = 'workding direcotry (default: /scratch.global/zhoux379)'
    )
    args = parser.parse_args()
    dirw = args.dirw

    print("working in %s" % dirw)
    os.chdir(dirw)
    ary = os.listdir(dirw)
    for item in ary:
        print(item)
        if item == 'temp' or item.startswith("sna"):
            continue
        cmd = "cp -rf %s %s.bak" % (item, item)
        print("  " + cmd)
        os.system(cmd)
        cmd = "rm -rf %s" % item
        print("  " + cmd)
        os.system(cmd)
        cmd = "mv %s.bak %s" % (item, item)
        print("  " + cmd)
        os.system(cmd)
