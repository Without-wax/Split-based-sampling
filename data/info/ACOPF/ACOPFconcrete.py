# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:22:52 2017

@author: jlc516
"""
from __future__ import division
from pyomo.environ import *

import numpy as np

def createACOPFModel(net):
    model = ConcreteModel()
    
    ######################################################################################
    #CREATE SETS   
    model.B = Set(initialize = net["bus"].index.tolist(), doc='Buses')
    model.G = Set(initialize = net["gen"].index.tolist(), doc='Generators')
    model.LO = Set(initialize = net["load"].index.tolist(), doc='Loads')
    model.LI = Set(initialize = net["line"].index.tolist(), doc='Lines')  
    
    
    ######################################################################################
    #CREATE PARAMETERS
    Load_parameters = list(net["load"].axes[1])
    for item in Load_parameters:
        setattr(model, "lo_"+item, Param(model.LO, initialize = dict(net["load"][item])))
    
    Gen_parameters = list(net["gen"].axes[1])
    for item in Gen_parameters:
        setattr(model, "ge_"+item, Param(model.G, initialize = dict(net["gen"][item])))

    Line_parameters = list(net["line"].axes[1])
    for item in Line_parameters:
        setattr(model, "li_"+item, Param(model.LI, initialize = dict(net["line"][item])))

    Bus_parameters = list(net["bus"].axes[1])
    for item in Bus_parameters:
        setattr(model, "bu_"+item, Param(model.B, initialize = dict(net["bus"][item])))

    #scalars
    model.base = Param(initialize=net["par"]["base"][0])
    model.ref = Param(initialize=net["par"]["refnode"][0])    
    model.vadegree = Param(initialize=net["par"]["va_degree"][0])
    
    
    ######################################################################################
    #CREATE VARIABLES

    #generation
    def fga(model, i): return (model.ge_min_p[i], model.ge_max_p[i])
    def iga(model,i): return model.ge_ini_p[i]
    model.ga = Var(model.G, initialize = iga, bounds = fga)
    
    def fgr(model, i): return (model.ge_min_q[i], model.ge_max_q[i])
    def igr(model,i): return model.ge_ini_q[i]
    model.gr = Var(model.G, initialize = igr, bounds = fgr)
    
    #voltages for all buses
    def vt(model, i):
        for k in model.G:
            if model.ge_bus[k] == model.bu_name[i]: 
                #model.volt[i].fixed = True
                return model.ge_vg[k]
        return 1
    
    def fv(model,i): return (model.bu_v_min[i], model.bu_v_max[i])
    model.volt = Var(model.B, initialize = vt, bounds = fv)
      
    
    #phi for all buses
    def ph(model, i):
        if model.bu_name[i] == net["par"]["refnode"][0]:
            model.phi[i].fixed = True 
            return net["par"]["va_degree"][0]
        else: return 0
    model.phi = Var(model.B, initialize = ph)
             
    #line flow
    #def ff(model, i): return (- model.li_max_f[i], model.li_max_f[i])
    model.fas = Var(model.LI)
    model.far = Var(model.LI)
    model.frs = Var(model.LI)
    model.frr = Var(model.LI)
    
    
    ######################################################################################
    #CREATE Constraints
    def node_balance_active(model, b):
        return - sum(model.ga[k] for k in model.G if model.ge_bus[k] == model.bu_name[b]) \
            + sum(model.lo_p[k] for k in model.LO if model.lo_bus[k] == model.bu_name[b]) \
            + sum(model.fas[k] for k in model.LI if model.li_from_bus[k] == model.bu_name[b]) \
            + sum(model.far[k] for k in model.LI if model.li_to_bus[k] == model.bu_name[b]) \
            == 0
    model.NodeBalanceActive = Constraint(model.B, rule=node_balance_active)
    
    
    def node_balance_reactive(model, b):
        return - sum(model.gr[k] for k in model.G if model.ge_bus[k] == model.bu_name[b]) \
            + sum(model.lo_q[k] for k in model.LO if model.lo_bus[k] == model.bu_name[b]) \
            + sum(model.frs[k] for k in model.LI if model.li_from_bus[k] == model.bu_name[b]) \
            + sum(model.frr[k] for k in model.LI if model.li_to_bus[k] == model.bu_name[b]) \
            == 0
    model.NodeBalanceReactive = Constraint(model.B, rule=node_balance_reactive)
    
    
    def line_active_send(model, l):
        return model.fas[l] == model.base * ( model.li_g[l] * sum(model.volt[k] * model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k] ) \
                       - sum(model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) * sum(model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) * ( \
                       + model.li_g[l] * cos( sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) ) \
                       + model.li_b[l] * sin( sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) ) ) )
    model.LinePhaseActiveSend = Constraint(model.LI, rule=line_active_send)
    
    def line_active_receive(model, l):
        return model.far[l] == model.base * ( model.li_g[l] * sum(model.volt[k] * model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k] ) \
                       - sum(model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) * sum(model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) * ( \
                       + model.li_g[l] * cos( sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) ) \
                       + model.li_b[l] * sin( sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) ) ) )
    model.LinePhaseActiveReceive = Constraint(model.LI, rule=line_active_receive)
    
    def line_reactive_send(model, l):
        return model.frs[l] == model.base * ( - (model.li_b[l] + model.li_b_cap[l]/2)  * sum(model.volt[k] * model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k] ) \
                       + sum(model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) * sum(model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) * ( \
                       - model.li_g[l] * sin( sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) ) \
                       + model.li_b[l] * cos( sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) ) ) )
    model.LinePhaseReactiveSend = Constraint(model.LI, rule=line_reactive_send)
    
    def line_reactive_receive(model, l):
        return model.frr[l] == model.base * ( - (model.li_b[l] + model.li_b_cap[l]/2)   * sum(model.volt[k] * model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k] ) \
                       + sum(model.volt[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) * sum(model.volt[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) * ( \
                       - model.li_g[l] * sin( sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) ) \
                       + model.li_b[l] * cos( sum(model.phi[k] for k in model.B if model.li_to_bus[l] == model.bu_name[k]) - sum(model.phi[k] for k in model.B if model.li_from_bus[l] == model.bu_name[k]) ) ) )
    model.LinePhaseReactiveReceive = Constraint(model.LI, rule=line_reactive_receive)
    
    
    def line_flow_send(model, l):
        return (0 , model.fas[l]*model.fas[l] + model.frs[l]*model.frs[l], model.li_max_f[l]*model.li_max_f[l])
    #model.LineFlowSend = Constraint(model.LI, rule=line_flow_send)
    
    def line_flow_receive(model, l):
        return (0 , model.far[l]*model.far[l] + model.frr[l]*model.frr[l], model.li_max_f[l]*model.li_max_f[l])
    #model.LineFlowReceive = Constraint(model.LI, rule=line_flow_receive)
 
    ######################################################################################
    #CREATE Constraints
    def objective_mincost(model): return sum(model.ge_cost_a[k] + model.ga[k] * model.ge_cost_b[k] + model.ga[k] * model.ga[k] * model.ge_cost_c[k] for k in model.G)
    model.objective = Objective(rule=objective_mincost, sense=minimize)
    
    return model