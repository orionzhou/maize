# -*- coding: utf-8 -*-
import os, sys, getopt;
import datetime;
import my;

fNameLst = ["hap{0}.txt", "hapinfo{0}.txt", "haploview{0}.txt", "rsq{0}.out"];
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    newFileName = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            newFileName = a
        else:
            assert False
    if newFileName == None:
        now = datetime.datetime.now();
        newFileName = now.strftime("%Y-%m-%d_%H:%M");
    for fName in fNameLst:
        fIn = os.path.join(my.DIR_stat, fName.format(""));
        fOut = os.path.join(my.DIR_site, "LD", fName.format("_"+newFileName));
        os.system("cp {0} {1}".format(fIn, fOut));

def usage():
    print "Usage: {0} -o <output_file_name>".format(sys.argv[0]);
    
if __name__ == "__main__":
    main()
