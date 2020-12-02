# Monte Calro integral
import random as rd

def f(x):
    return 4/(1+x*x)
	
N = 1000000
I = 0
for n in range(N):
    I += f(rd.random())

print(I/N)
