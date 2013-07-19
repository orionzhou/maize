# -*- coding: utf-8 -*-
import sys
import os
import re

DIR_Src = os.environ["src"]
DIR_Mt3 = os.path.join(os.environ["m"], "mt3")
DIR_Work = os.environ["work"]
DIR_Data = os.environ["data"]
DIR_Genome = os.path.join(DIR_Data, "genome")
DIR_In = os.path.join(DIR_Data, "in")
DIR_Out = os.path.join(DIR_Data, "out")
DIR_Misc1 = os.path.join(DIR_Data, "misc1")
DIR_Misc2 = os.path.join(DIR_Data, "misc2")
DIR_Misc3 = os.path.join(DIR_Data, "misc3")
DIR_Tmp = "/project/scratch/zhoup/tmp"

if __name__ == "__main__":   
    print DIR_Work, DIR_Data
