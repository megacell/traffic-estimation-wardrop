'''
Created on Apr 23, 2014

@author: jeromethai
'''

import path_solver as path
import test_ue_solver as testue
import missing as mis
import numpy as np
import matplotlib.pyplot as plt


def main():
    grid, linkflows, unusedpaths = testue.main()
    errors = mis.avg_error(grid, [0, 1, 2, 3, 4, 5, 6, 7], 50)
    x = np.arange(8)
    xlabel = []
    for i in x: xlabel.append('%i%%' % ((100*i)/8))
    fig, ax = plt.subplots()
    plt.bar(x, errors)
    plt.xticks( x + 0.5,  xlabel )
    plt.xlabel('Percentage of missing values')
    plt.ylabel('Error in link flows')
    plt.title(r'Relative $\ell$-2 error for link flows')
    plt.show()

    
if __name__ == '__main__':
    main()