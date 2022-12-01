#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 10:21:23 2022
 decided to use lorenz attractor because it is so well documented
on scale from one to ten this is ugly
@author: CJC
"""
# Quick Attractor Simulation For Random
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

def change_coeff(coeff, sigma, rho, beta):
    coeff.update({'sigma':sigma})
    coeff.update({'rho':rho})
    coeff.update({'beta':beta})

def change_states(states, x, y, z):
    states.update({'x':x})
    states.update({'y':y})
    states.update({'z':z})

def dx_dt(sigma, x, y):
    dxdt = sigma*(y-x)
    return dxdt

def dy_dt(rho, x, y, z):
    dydt = rho*x - y - x*z
    return dydt

def dz_dt(beta, x, y, z):
    dzdt = x*y -beta*z
    return dzdt
                        # coefficients
coeff = {'sigma': 10.,
         'rho': 28.000,
         'beta': 2.667}

                        # states stored in dictionary - altered over time
states = {'x': 0.,
          'y': 1.,
          'z': 1.05,
          't': 0}

                        # initial states, just in case
istates = {'x': 0.,
          'y': 1.,
          'z': 1.05,
          't': 0}

                        # states in dictionary altered each instance
                        # all parameters here can be swapped out as long as they are provided as dictionaries
                        
def update_states(dt = 0.01, states = states, coeff = coeff):
    dxdt = dx_dt(coeff.get('sigma'), states.get('x'), states.get('y'))
    dydt = dy_dt(coeff.get('rho'), states.get('x'), states.get('y'), states.get('z'))
    dzdt = dz_dt(coeff.get('beta'), states.get('x'), states.get('y'), states.get('z'))
    t = (states.get('t') + 1)
    x_state = states.get('x') + (dxdt*dt)
    y_state = states.get('y') + (dydt*dt)
    z_state = states.get('z') + (dzdt*dt)
    states.update({'x':x_state})
    states.update({'y':y_state})
    states.update({'z':z_state})
    states.update({'t':t})
    return [x_state, y_state, z_state, t] # returns a state
                        # 
def run_slice(time, states, coeff, dt = 0.01):
    outputs = list()
    for t in range(0,time):
        state = update_states(dt, states, coeff)
        outputs.append(state)
    return outputs

def wtimeseries(slice_data):
    with open('tempdata.csv', 'w') as tempwrite:
        for line in slice_data:
            for subline in line:
                tempwrite.write(f'{subline},')
            tempwrite.write(f'\n')

def show_slice(slice_data):
    pio.renderers.default='browser'
    x=np.array([sub[0] for sub in slice_data])
    y=np.array([sub[1] for sub in slice_data])
    z=np.array([sub[2] for sub in slice_data])
    fig = go.Figure(data=[go.Scatter3d(
        x=x,
        y=y,
        z=z
        )])

    # tight layout
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.show()










