import numpy as np
import numpy.linalg as alg
import numpy.random as r
import scipy.linalg as alg
import scipy.sparse as spa
from time import time

def Torus_Cov( n, iter ):
# Ce script calcule une sequence de temp couverture en une torus annulus 
# T_N par MCMC

    # Parameters for random walk
    X_t     = 0
    Y_t     = 0
    T_cov   = np.zeros(iter)

    for i in range(iter):
        L_t     = np.zeros([n, n])            # Local time matrix
        n_V     = n * n                       # of not-v node
        t       = 0
        while n_V != 0:
            if(L_t[X_t, Y_t] == 0):
                n_V = n_V - 1
                L_t[X_t, Y_t] = 1
            t   = t + 1
            rnd = r.rand(1)
            if rnd > 0.75:
                X_t = np.mod(X_t + 1, n)
            elif rnd > 0.5:
                X_t = np.mod(X_t - 1, n)
            elif rnd > 0.25:
                Y_t = np.mod(Y_t + 1, n)
            else:
                Y_t = np.mod(Y_t - 1, n)
        print(t / np.power(n * np.log(n), 2) * np.pi / 4)
        T_cov[i]    = t
        
    return T_cov


start = time()
n = 50
iter = 100
T_cov = Torus_Cov(n, iter)
print(time() - start)
write = False
if write:
    with open('50.txt','a',encoding='utf-8') as f:#使用with open()新建对象f
        for i in range(iter):
            f.write(str(T_cov[i]))
            f.write('\n')