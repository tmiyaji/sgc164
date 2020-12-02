# compute the area of the unit circle by Monte Calro
import random as rd
		
N = 1000000
S = 4
M = 0
for n in range(N):
    x = -1 + 2 * rd.random()
    y = -1 + 2 * rd.random()
    M += 1 if x*x+y*y < 1 else 0

print(S*M/N)
