#!/usr/bin/env python
# coding: utf-8


# !pip install pandapower==2.6.0
# print(pp.__version__)
import numpy as np
import pandas as pd
import pandapower as pp
import pandapower.networks as pn
import pandapower.converter as pc
from pypower.api import makeYbus
from pyomo.opt import SolverFactory
from pyomo.environ import *
from pyomo.environ import value
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import time
import cmath
# from pandapower.results import calculate_voltage_dependent_losses

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

def _maxgap(points, input_space=None):
    # points has rows = samples, columns = variables

    # make a copy before sorting the columns
    points = points.copy()
    if input_space is not None:
        points = np.vstack((input_space.minimum, points, input_space.maximum))
    points.sort(0)

    gaps = points[1:,:] - points[:-1,:]
    width = gaps.max(0)
    loc = gaps.argmax(0)
    left = points[loc, np.arange(points.shape[1])]
    relative = width / (points[-1,:] - points[0,:])
    target = left + width/2

    return relative, target, width


def _generate_Ybus(net, baseMVA, path = 'data/case14.xlsx'):
# pp.runpp(net)
# ppc = pp.converter.to_ppc(net)
# Ybus = net["_ppc"]["internal"]["Ybus"]
# Using pandas
# Ybus, Yf, Yt = makeYbus(baseMVA, ppc['bus'], ppc['branch'])
# bus_len = len(Ybus)

    xls = pd.ExcelFile(path)
    branch = pd.read_excel(xls, 'branch')
    gen = pd.read_excel(xls, 'gen')
    bus = pd.read_excel(xls, 'bus')
    bus_len = bus.shape[0]; gen_len = len(gen['bus'].ravel())

    Ybus = np.matrix(np.zeros([bus_len,bus_len]),dtype=complex)
    shunt = np.zeros(bus_len, dtype=complex)
    x = [pd.DataFrame() for i in range(bus_len)]
    for i in range(bus_len):  
        x[i] = branch.loc[branch['fbus'] == i]
        fro = [j for j in x[i]['tbus']]
        fro = pd.DataFrame(fro)
        for j in range(len(fro)):
            Ybus[i-1,fro.loc[j][0]-1] =  -1/complex(pd.DataFrame(x[i]['r']).iloc[j][0],pd.DataFrame(x[i]['x']).iloc[j][0])
            Ybus[fro.loc[j][0]-1,i-1] =  -1/complex(pd.DataFrame(x[i]['r']).iloc[j][0],pd.DataFrame(x[i]['x']).iloc[j][0])
            shunt[j] = complex(0,sum(x[i]['b'])/2)
    for i in range(bus_len):
        Ybus[i,i] = np.sum(Ybus[i,:])*(-1) + complex(0,sum(branch.loc[(branch['fbus']==i+1) | (branch['tbus']==i+1)]['b'])/2)
        
    return Ybus, bus_len, gen_len



def _get_network_properties(net, baseMVA):
    #Get network elements
    branch = net.line
    gen = net.gen
    bus = net.bus
    load = net.load
    cost = net.poly_cost

    #Get all generators (+ external grid)
    comb_gen = gen[list(net.ext_grid.columns.drop('va_degree'))]
    comb_gen = pd.concat([comb_gen, net.ext_grid]).sort_values(by=['bus'])
    gen_list = list(comb_gen.bus)
    bus_len = len(bus)
    
    #Calculate Max MVA for branches
    max_mva = []
    for i in range(len(branch)):
        max_mva.append(net.line.max_i_ka[i]*np.sqrt(3)*net.bus.iloc[int(net.line.from_bus[i])]['vn_kv']) 
    branch['rateA'] = max_mva

    # initialisation for input space
    Pd = np.zeros(bus_len); Qd = np.zeros(bus_len)
    Vmax = np.zeros(bus_len); Vmin = np.zeros(bus_len)
    Pmax = np.zeros(bus_len); Pmin = np.zeros(bus_len)
    Qmax = np.zeros(bus_len); Qmin = np.zeros(bus_len)
#     cost_lin = np.zeros(bus_len); cost_quad = np.zeros(bus_len)

    for i in range(bus_len):
        Vmax[i] = bus['max_vm_pu'][i]
        Vmin[i] = bus['min_vm_pu'][i]
    for i,j in enumerate(load.bus):
        Pd[j] = load.p_mw.iloc[i]/baseMVA
        Qd[j] = load.q_mvar.iloc[i]/baseMVA

    for i,j in enumerate(comb_gen.bus):
        Pmax[j] = comb_gen.max_p_mw.iloc[i]/baseMVA
        Pmin[j] = comb_gen.min_p_mw.iloc[i]/baseMVA
        Qmax[j] = comb_gen.max_q_mvar.iloc[i]/baseMVA
        Qmin[j] = comb_gen.min_q_mvar.iloc[i]/baseMVA
#         cost_lin[j] = cost.cp1_eur_per_mw.iloc[i]*baseMVA
#         cost_quad[j] = cost.cp2_eur_per_mw2.iloc[i]*baseMVA**2


    cols = ['minimum','maximum']
    p = comb_gen[['min_p_mw','max_p_mw']]/baseMVA
    p.columns = cols
    p = p.T
    p.columns = list(comb_gen.reset_index().bus)
    # p.columns = ['P%s' % (i+1) for i in range(len(comb_gen))]
    p = p.T 

    q = comb_gen[['min_q_mvar','max_q_mvar']]/baseMVA
    q.columns = cols
    q = q.T
    q.columns = list(comb_gen.reset_index().bus)
    # q.columns = ['Q%s' % (i+1) for i in range(len(comb_gen))]
    q = q.T 

    # input_space = pd.concat([p,q])
    input_space = p

    return input_space, Pmax, Pmin, Qmax, Qmin, Vmax, Vmin, Pd, Qd, branch, gen_list



def _generate_sample(
        primary_var, primary_lb, primary_ub,
        secondary_vars, secondary_targets, secondary_weights,
        gen_len, bus_len, Vmax, Vmin, Pmax, Pmin, Qmax, Qmin,
        input_space, Ybus, branch, Pd, Qd, gen_list):

    
    #Create a model and declare the variables
    model = ConcreteModel()
    model.IDX1 = range(gen_len)
    model.IDX2 = range(bus_len)
    model.Qg = Var(model.IDX2)
    model.t = Var(model.IDX2, initialize=[0 for i in model.IDX2])
    model.v = Var(model.IDX2, initialize=[1 for i in model.IDX2])
    model.vars = Var(model.IDX2)

    model.c = ConstraintList()

    #Bus voltage angle initial condition
    model.c.add(model.t[0] == 0)
#     Bounds on voltage, reactive and acitve power generation
    for i in range(bus_len):
        model.c.add(expr = model.v[i] <= Vmax[i])
        model.c.add(expr = model.v[i] >= Vmin[i])
        model.c.add(expr = model.Qg[i] <= Qmax[i])
        model.c.add(expr = model.Qg[i] >= Qmin[i])
        if i in gen_list:
            model.Qg[i].fixed = False
            model.c.add(expr = model.vars[i] <= input_space.loc[i].maximum)
            model.c.add(expr = model.vars[i] >= input_space.loc[i].minimum)
        else:        
            model.Qg[i].fix(0)
            model.vars[i].fix(0)
            
#     set primary variable
    model.vars[primary_var].setlb(primary_lb)
    model.vars[primary_var].setub(primary_ub)

    if secondary_vars is not None:
        quad_exp = 0
        for i, sec in enumerate(secondary_vars):
            quad_exp += secondary_weights[i] * (model.vars[input_space.index[sec]]-secondary_targets[i])**2

        #Nodal equations
#         count = 0
        for i in range(bus_len):
            model.c.add(expr = sum([model.v[i]*model.v[j]*cmath.polar(Ybus[i,j])[0]*cos(model.t[i]-model.t[j]-cmath.polar(Ybus[i,j])[1]) for j in range(bus_len)]) - model.vars[i] +  Pd[i] == 0)
            model.c.add(expr = sum([model.v[i]*model.v[j]*cmath.polar(Ybus[i,j])[0]*sin(model.t[i]-model.t[j]-cmath.polar(Ybus[i,j])[1]) for j in range(bus_len)]) - model.Qg[i] +  Qd[i] == 0)

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

        model.cost = Objective(expr = quad_exp)
    else:
        model.cost = Objective(expr = 1.0)
    
    solution = SolverFactory('ipopt').solve(model)

    if (solution.solver.status == SolverStatus.ok) and (solution.solver.termination_condition == TerminationCondition.optimal):
        # Do something when the solution in optimal and feasible
        model.vars[primary_var].setlb(input_space.iloc[primary_var]['minimum'])
        model.vars[primary_var].setub(input_space.iloc[primary_var]['maximum'])

        result = np.array([model.vars[i].value for i in list(input_space.index)])
        result1 = np.array([model.Qg[i].value for i in list(input_space.index)])
        return result, result1

    else:
        model.vars[primary_var].setlb(input_space.iloc[primary_var]['minimum'])
        model.vars[primary_var].setub(input_space.iloc[primary_var]['maximum'])
        return None, None


def GS(net, _path = 'data/case14.xlsx', n=100, secondary_frac=0.5, start_time = time.time(), report_interval=0.1): 
    #initialisations
    k = 0
    infeasible_count = 0
    true_count = 0
    primary_var = -1
    try_ = 0
    max_tries = None; sample_upper_bound = 10**6; min_range=1e-5; primary_tol=0.001; enforce_range=True,
    report = lambda s: print(s)

    #Calculate Ybus
    Ybus, bus_len, gen_len = _generate_Ybus(net, net.sn_mva, path = _path)
    #Compute bounds of the input space
    input_space, Pmax, Pmin, Qmax, Qmin, Vmax, Vmin, Pd, Qd, branch, gen_list = _get_network_properties(net, net.sn_mva)

#     for secondary variables
    if secondary_frac >= 1.0:
        n_secondary = secondary_frac
    else:
        n_secondary = np.floor(secondary_frac * len(input_space)).astype(int)
#     get indices and weights of variables
    idxs = (input_space.maximum - input_space.minimum >= min_range).to_numpy().nonzero()[0]
    weights = (1/(input_space.maximum - input_space.minimum)**2).to_numpy()
#     report output
    report("Targeting {}/{} unblocked primary variables.".format(len(idxs), len(input_space)))
    report("Targeting {} secondary variables.".format(n_secondary))
    report_header, report_format = _make_report_header(n)
    report("\n" + report_header)
    if report_interval < 1.0:
        report_interval = np.floor(report_interval * n).astype(int)

    #arrays of all generated samples, all feasible samples and corresponding reactive set-points
    samples = np.zeros((sample_upper_bound, len(input_space)))
    i_samples = np.zeros((n, len(input_space)))
    reactive_samples = np.zeros((sample_upper_bound, len(input_space)))
    #index of primary variables
    var_p = np.array([]); 
    #index of feasible samples
    sample_id = np.array([])

    while True:
        if max_tries is not None and try_ >= max_tries:
            break
        try_ += 1

        relative, target, width = _maxgap(samples[0:k,idxs], input_space.iloc[idxs,:])
    #     Primary variable as the index with the largest gap
        primary_var = np.argmax(relative)

        #Selected primary variable
        var_p = np.append(var_p,idxs[primary_var])
#         identify primary and secondary variables, targets
        primary_target = target[primary_var]
        primary_lb = primary_target - primary_tol*width[primary_var]
        primary_ub = primary_target + primary_tol*width[primary_var]
        secondary_vars = np.random.choice(len(idxs), n_secondary, replace=False)
        secondary_targets = target[secondary_vars]
        secondary_weights = weights[idxs[secondary_vars]]

        #generate new sample
        new_sample, reactive =  _generate_sample(
                    idxs[primary_var], primary_lb, primary_ub,
                    idxs[secondary_vars], secondary_targets, secondary_weights,
                    gen_len, bus_len, Vmax, Vmin, Pmax, Pmin, Qmax, Qmin,
                    input_space, Ybus, branch, Pd, Qd, gen_list)

        if new_sample is not None:
            if enforce_range:
                new_sample[new_sample > input_space.maximum] = input_space.maximum[new_sample > input_space.maximum]
                new_sample[new_sample < input_space.minimum] = input_space.minimum[new_sample < input_space.minimum]
            # Store generated sample
            samples[k,:] = new_sample
            i_samples[true_count,:] = new_sample
            # Store reactive power set-points Qg
            reactive_samples[k,:] = reactive

            sample_id = np.append(sample_id,k)
            k += 1
            true_count += 1


            #Reporting true samples
            relative, target, width = _maxgap(i_samples[0:true_count,idxs], input_space.iloc[idxs,:])
            if k % report_interval == 0:
                elapsed = time.time() - start_time
                remaining = ((n - true_count) / n )*100
                report(report_format.format(
                        i=true_count, n=n, ela=elapsed, rem=remaining,
                        inf=infeasible_count))
        else:
    #         Store the infeasible sample for primary variable without altering the gaps in other variables
            samples[k,primary_var] = target[primary_var]
            w = idxs.copy()
            w = w[w != idxs[primary_var]]
            for i in range(len(w)):
                samples[k,w[i]] = input_space['minimum'].iloc[w[i]]
            k += 1
            infeasible_count += 1         

        if true_count >= n: break

    if true_count < n:
        # max_tries reached; return fewer than n samples
        i_samples = i_samples[:true_count,:]  

    if k < sample_upper_bound:
        # max_tries reached; return fewer than n samples
        samples = samples[:k,:]
        reactive_samples = reactive_samples[:k,:]
    
    all_samples = pd.DataFrame(data=samples,columns=input_space.maximum.index)
    feasible_samples = pd.DataFrame(data=i_samples,columns=input_space.maximum.index)
    reactive_powers = pd.DataFrame(data=reactive_samples, columns=input_space.maximum.index)
    
    return all_samples, feasible_samples, reactive_powers

