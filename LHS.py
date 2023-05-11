#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
from smt.sampling_methods import LHS

import pandapower as pp
import pandapower.networks as pn
import pandapower.converter as pc
from pypower.api import makeYbus
from pyomo.opt import SolverFactory
from pyomo.environ import *
from pyomo.environ import value
from GAPSPLIT import  _generate_Ybus, _get_network_properties
import time
import cmath

def _make_report_header(maxN):
    """Return the header and format string for reporting coverage."""
    nw = len(str(maxN))
    frac_width = 2*nw + 1  # width of 300/1000
    frac_header = 'Sample'
    frac_format = '{i:' + str(nw) + 'd}/{n:' + str(nw) + 'd}'
    if frac_width < len(frac_header):
        pad = ''.join([' ' for _ in range(len(frac_header) - frac_width)])
        frac_format = pad + frac_format
    elif len(frac_header) < frac_width:
        pad = ''.join([' ' for _ in range(frac_width - len(frac_header))])
        frac_header = pad + frac_header
    hdr = frac_header + " Elapsed (s)     Remaining(%)   Infeasible"
    fmt = frac_format + " {ela:9.2f}     {rem:9.2f}   {inf:10d}"
    return hdr, fmt
    
def _generate_sample_RS(pg, gen_len, bus_len, Vmax, Vmin, Pmax, Pmin, Qmax, Qmin,
        input_space, Ybus, branch, Pd, Qd, gen_list):#qg
        #Create a model
    model = ConcreteModel()

    model.IDX1 = range(gen_len)
    model.IDX2 = range(bus_len)
    model.Pg = Var(model.IDX2)
    model.Qg = Var(model.IDX2)
    model.t = Var(model.IDX2, initialize=[0 for i in model.IDX2])
    model.v = Var(model.IDX2, initialize=[1 for i in model.IDX2])
    
    #Create slack variables
    model.d1 = Var(model.IDX2)
#     model.d2 = Var(model.IDX2) 
     
    # Fix data from generated samples
    for i in model.IDX2:
        model.Pg[i].fix(0) 
        model.Qg[i].fix(0)
        model.d1[i].fix(0)
        if i in gen_list:
            model.Qg[i].fixed = False
            model.d1[i].fixed = False
            model.Pg[i].fix(float(pg.loc[i]))
        
    # declare constraints
    model.c = ConstraintList()
    
    for i in model.IDX2:
        model.c.add(expr = model.v[i] <= Vmax[i])
        model.c.add(expr = model.v[i] >= Vmin[i])
        model.c.add(expr = model.Qg[i] <= Qmax[i])
        model.c.add(expr = model.Qg[i] >= Qmin[i])
    
        if i in gen_list:
            model.c.add(expr = (model.Pg[i] + model.d1[i]) <= Pmax[i])
            model.c.add(expr = (model.Pg[i] + model.d1[i]) >= Pmin[i])

    #Bus voltage angle initial condition
    model.c.add(model.t[0] == 0)

    #Nodal equations
    for i in range(bus_len):
        model.c.add(expr = sum([model.v[i]*model.v[j]*cmath.polar(Ybus[i,j])[0]*cos(model.t[i]-model.t[j]-cmath.polar(Ybus[i,j])[1]) for j in range(bus_len)]) - model.Pg[i] + Pd[i] == model.d1[i])
        model.c.add(expr = sum([model.v[i]*model.v[j]*cmath.polar(Ybus[i,j])[0]*sin(model.t[i]-model.t[j]-cmath.polar(Ybus[i,j])[1]) for j in range(bus_len)]) - model.Qg[i] + Qd[i] == 0)#model.d2[i]

    #Line flow constraints
    for i in range(len(branch)): 
        x = branch.from_bus[i]
        y = branch.to_bus[i]
        val = branch.rateA[i]/100
        if(val == 0):
            val = 100
    #    With phasor
        Pxy = (model.v[x]**2)*cmath.polar(Ybus[x,y])[0]*cos(cmath.polar(Ybus[x,y])[1]) - model.v[x]*model.v[y]*cmath.polar(Ybus[x,y])[0]*cos(model.t[x]-model.t[y]-cmath.polar(Ybus[x,y])[1])
        Qxy =-(model.v[x]**2)*cmath.polar(Ybus[x,y])[0]*sin(cmath.polar(Ybus[x,y])[1]) - model.v[x]*model.v[y]*cmath.polar(Ybus[x,y])[0]*sin(model.t[x]-model.t[y]-cmath.polar(Ybus[x,y])[1])     

        Pyx = (model.v[y]**2)*cmath.polar(Ybus[x,y])[0]*cos(cmath.polar(Ybus[x,y])[1]) - model.v[x]*model.v[y]*cmath.polar(Ybus[x,y])[0]*cos(model.t[y]-model.t[x]-cmath.polar(Ybus[x,y])[1])
        Qyx =-(model.v[y]**2)*cmath.polar(Ybus[x,y])[0]*sin(cmath.polar(Ybus[x,y])[1]) - model.v[x]*model.v[y]*cmath.polar(Ybus[x,y])[0]*sin(model.t[y]-model.t[x]-cmath.polar(Ybus[x,y])[1])     

        if (value(Pxy) != 0) & (value(Pyx) != 0) :
            model.c.add(expr = (Pxy**2 + Qxy**2)<= val**2)
            model.c.add(expr = (Pyx**2 + Qyx**2)<= val**2)


 #Declare objective function
    model.cost = Objective( 
        expr = sum([model.d1[i]**2 for i in model.IDX2]) #+ sum([model.d2[i]**2 for i in range(gen_len)]) 
    ) 

    # solve
    solution = SolverFactory('ipopt').solve(model)
   
    
    if (solution.solver.status == SolverStatus.ok) and (solution.solver.termination_condition == TerminationCondition.optimal):
        if(model.cost()  >= 0.01):
            result = np.array([float(model.Pg[i].value + model.d1[i].value) for i in gen_list]
#                               +[(model.Qg[i].value + model.d2[i].value) for i in gen_list]
                             )
#             result = np.array([1 for i in model.IDX1])
#             print("Slight mod")
            return result
        else:
            result = np.array([float(model.Pg[i].value) for i in gen_list]
#                               +[(model.Qg[i].value) for i in gen_list]
                             )
            return result
        
    else:
        return np.array([1 for i in model.IDX1])

    
def RS(net, _path = 'data/case14.xlsx', n=100, start_time = time.time(), report_interval=0.1):  

    #Calculate Ybus
    Ybus, bus_len, gen_len = _generate_Ybus(net, net.sn_mva, path = _path)
    #Compute bounds of the input space
    input_space, Pmax, Pmin, Qmax, Qmin, Vmax, Vmin, Pd, Qd, branch, gen_list = _get_network_properties(net, net.sn_mva)
    sample_upper_bound = 10**6;
    pg = np.array(input_space)
    sampling_pg = LHS(xlimits=pg)
    gen_active = sampling_pg(sample_upper_bound)
#     qg = np.array(q)
#     sampling_qg = LHS(xlimits=qg)
#     gen_reactive = sampling_qg(10*n)
    report = lambda s: print(s)
    report_header, report_format = _make_report_header(n)
    report("\n" + report_header)
    if report_interval < 1.0:
        report_interval = np.floor(report_interval * n).astype(int)
        
    k = 0
    samples = np.zeros(gen_len)
    infeasible_samples = np.array([])
    true_count = 0;
    infeasible_count = 0;
    while(1):
        if (len(samples)-len(infeasible_samples)-1 == n):
            break
        sample = pd.DataFrame(gen_active[k,:])
        sample.index = input_space.index
        new_sample = _generate_sample_RS(sample, gen_len, bus_len, Vmax, Vmin, Pmax, Pmin, Qmax, Qmin,
        input_space, Ybus, branch, Pd, Qd, gen_list)#, gen_reactive[k,:]
        if new_sample is not None:
            if(np.all(new_sample == 1)):
                infeasible_samples = np.append(infeasible_samples,k)
                infeasible_count += 1
            else:
                true_count += 1
            samples = np.vstack((samples,new_sample))
            k += 1            
        if k % report_interval == 0:
            elapsed = time.time() - start_time
            remaining = ((n - true_count) / n )*100
            report(report_format.format(
                    i=k, n=n, ela=elapsed, rem=remaining,
                    inf=infeasible_count))
            
    df_samples = pd.DataFrame(data=samples[1:], columns=input_space.maximum.index)
    df_samples = df_samples.drop(infeasible_samples)
#     df_samples.columns = ['P%s' %(i+1) for i in range(gen_len)]#+['Q%s' %(i+1) for i in range(gen_len)]
    return df_samples, pd.DataFrame(gen_active).drop(infeasible_samples).iloc[0:n]

