#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 19:05:25 2026

@author: shiz16
"""
import numpy as np
from scipy.optimize import minimize, LinearConstraint, Bounds

def EIOT_PRED(S_E, NUM_NC, dm, c_A_bounds, beq):

    
    # Compute H and f matrices
    H = S_E.T @ S_E
    f = -S_E.T @ dm
    
    A = np.array([]).reshape(0, S_E.shape[1])  # Empty constraint matrix
    b = np.array([])
    
    # Set bounds and equality constraints
    if len(c_A_bounds) != 0:
        # Constrained case
        Aeq = np.hstack([np.ones((1, S_E.shape[1] - NUM_NC)), 
                          np.zeros((1, NUM_NC))])
        lb = np.hstack([np.zeros(S_E.shape[1] - NUM_NC), 
                        np.ones(NUM_NC) * c_A_bounds[0]])
        ub = np.hstack([np.ones(S_E.shape[1] - NUM_NC), 
                        np.ones(NUM_NC) * c_A_bounds[1]])
    else:
        # Unconstrained case
        Aeq = np.ones((1, S_E.shape[1]))
        lb = np.zeros(S_E.shape[1])
        ub = np.ones(S_E.shape[1])
    
    # Define objective function: 0.5 * x^T @ H @ x + f @ x
    def objective(x):
        return 0.5 * x @ H @ x + f @ x
    
    # Define equality constraint: Aeq @ x = beq
    constraints = LinearConstraint(Aeq, beq, beq)
    bounds = Bounds(lb, ub)
    
    # Solve quadratic program
    result = minimize(
        objective,
        x0=np.zeros(S_E.shape[1]),
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={'ftol': 1e-9, 'maxiter': 1000}
    )
    
    c_E_hat = result.x
    fval = result.fun
    exitflag = result.success
    output = result
    lambda_vals = result
    
    # Compute results
    dm_hat = S_E @ c_E_hat
    Em = dm - dm_hat
    sse = np.sum((dm - dm_hat) ** 2)
    
    return c_E_hat, Em, sse, lambda_vals


"""
test code below
"""

import scipy.io as sp
data = sp.loadmat('Form_conversion_test.mat')

# Extract the necessary matrices
K = data['K_SG']  # K is the pure component spectra matrix (m x lambda)
X = data['X_SG']  # X is the mixture spectra matrix (n x lambda)


S_E = np.transpose(K)
NUM_NC = 0
dm = np.transpose(X[0,:])
c_A_bounds = [0, 1]
beq = 1.0

c_E_hat, Em, sse, lambda_vals = EIOT_PRED(S_E, NUM_NC, dm, c_A_bounds, beq)

print(c_E_hat)
