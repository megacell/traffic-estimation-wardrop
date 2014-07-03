'''
Created on Apr 23, 2014

@author: jeromethai
'''
import scipy.linalg as sla


def place_zeros(M, atol=1e-13):
    """Replace entries in M less than atol by 0.0"""
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if abs(M[i][j]) < atol: M[i][j] = 0.0
    return M


def find_basis(M):
    """Find the indices of the columns of M that form a basis or range(M)"""
    p,l,u = sla.lu(M)
    ind = [i for i in range(u.shape[0]) if u[i,i] != 0.0]
    if u[i,i] == 0:
        for j in range(i+1,u.shape[1]):
            if u[i,j] != 0.0: ind.append(j); break
    return ind


def bisection(F, f, left, right, tol=1e-8):
    """Use bisection to find x such that F(x)=f
    we suppose F is strictly increasing"""
    l, r = left, right
    while r-l>tol:
        if F((l+r)/2.0) < f:
            l = (l+r)/2.0
        else:
            r = (l+r)/2.0
    return (l+r)/2.0


if __name__ == '__main__':
    def F(x):
        return x**2
    print bisection(F, 9.0, 0.0, 5.0)
    