'''
Created on Jul 27, 2014

@author: jeromethai
'''

import numpy as np
import matplotlib.pyplot as plt
from cvxopt import matrix, spdiag

def true_delays(delaytype):
    """Display true delay functions"""
    coef = matrix([0.0, 0.0, 0.0, 0.15])
    degree = len(coef)
    a, b = 3.5, 3.0
    xdata = np.linspace(0.0, 2.5, num=100)
    if delaytype == 'Polynomial':
        true_vals = [1+(coef.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
        title = 'True polynomial delay function'
    if delaytype == 'Hyperbolic':
        true_vals = [1 - a/b + a/(b-x) for x in xdata]
        title = 'True hyperbolic delay function'
    plt.plot( xdata, true_vals)
    plt.xlabel('Link flow (1000 veh/h)')
    plt.ylabel('Delay')
    plt.title(title)
    plt.legend()
    plt.show()
 
 
def lin_delays(coef, current, type):
    """Display different linearizations"""
    degree = len(coef)
    xdata = np.linspace(0.0, 2.5, num=100)
    x0 = [1+(coef.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    #g = (coef.T * matrix(np.power(current,range(degree))))[0]
    g = (coef.T * spdiag(range(1,degree+1)) * matrix(np.power(current,range(degree))))[0]
    d = 1+(coef.T * matrix(np.power(current,range(1,degree+1))))[0]
    x1 = [x for x in xdata if d+g*(x-current) >= 1.0]
    y1 = [d+g*(x-current) for x in x1]
    plt.plot(x1, y1, 'r', label='tangent')
    plt.plot(xdata, x0, 'b', label='true')
    plt.plot(xdata, [d]*100, 'g', label='constant')
    plt.xlabel('Link flow (1000 veh/h)')
    plt.ylabel('Delay')
    plt.title('Linearized '+type+' delay function at x={}'.format(current))
    plt.legend()
    plt.show()
    
    xdata = np.linspace(0.0, 2.5, num=100)
    x0 = [1+(coef.T * matrix(np.power(x,range(1,degree+1))))[0] for x in xdata]
    a = (coef.T*matrix(np.power(0.85, range(degree))))[0]
    y2 = [1+a*x for x in xdata]
    plt.plot(xdata, y2, 'r', label='linearization')
    plt.plot(xdata, x0, 'b', label='true')
    plt.xlabel('Link flow (1000 veh/h)')
    plt.ylabel('Delay')
    plt.title('Linearized '+type+' delay function')
    plt.legend()
    plt.show()

if __name__ == '__main__':

    type = 'Hyperbolic'; coef = matrix([0.331, 0.217, 0.0104, 2.33e-7, 8.95e-8, 0.0107])
    #type = 'Polynomial'; coef = matrix([0.0, 0.0, 0.0, 0.15])
    #true_delays(type)
    lin_delays(coef, 2.0, type)