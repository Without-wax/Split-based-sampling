import numpy as np
import pandas as pd
import time
import pandapower.networks as pn
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from GAPSPLIT import _generate_sample, GS
from LHS import _generate_sample_RS, RS


#Load network (Change 'net' also '_path' to reflect the network of interest)
net = pn.case14(); 
baseMVA = 100
net.sn_mva = baseMVA

#define number of samples to generate
n = 10;

if __name__ == "__main__":
    #Generate OCs using split-based sampling (all samples, feasible samples, reactive samples)
    als, fs, rs = GS(net, _path = 'data/case14.xlsx', 
                    n=n, secondary_frac=0.5, start_time=time.time(),report_interval=0.1)


    #Generate OCs using LHS sampling (feasible rs, proposed rs)
    rs, proposed = RS(net, _path = 'data/case14.xlsx', 
                    n=n, start_time=time.time(), report_interval=0.1)

