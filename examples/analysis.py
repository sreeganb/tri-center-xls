#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import sys
import seaborn as sns
import glob
sys.path.append('/home/sree/git/PMI_analysis/pyext/src/')
from analysis_trajectories import *

if __name__ == "__main__":
    top_dir = sys.argv[1]
    analysis_dir = os.path.join(top_dir, "analysis/")
    # check if the analysis directory exists
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    # naming of the traj directory
    
