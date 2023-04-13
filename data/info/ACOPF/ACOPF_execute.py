#! /usr/bin/env python

# Python imports
from pyomo import *
from pyomo.opt.base import SolverFactory
import ACOPFconcrete as ACOPF
import pandas as pd

#read data from excel
net = pd.read_excel( open( "IEEE6AC.xlsx",'rb') , sheetname = None )

#create model
inst = ACOPF.createACOPFModel(net)

#create solver
solver = SolverFactory('ipopt')

#solver model
solver.solve(inst)