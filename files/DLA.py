import numpy as np
import random as rd
import matplotlib.pyplot as plt

class DLA:
    def __init__(self, N, L):
        self.seeds = [{'x':0, 'y':0}]
        self.walkers = [{'x':0, 'y':0, 'd':0} for k in range(N-1)]
        for n in range(N-1):
            s = rd.random() * 2 * np.pi
            self.walkers[n]['x'] = int(L * np.cos(s))
            self.walkers[n]['y'] = int(L * np.sin(s))
            self.walkers[n]['d'] = abs(self.walkers[n]['x']) + abs(self.walkers[n]['y'])

    def evolve(self):
        for walker in self.walkers:
            c = rd.randrange(0,4)
            if c == 0:
                walker['x'] += 1
            elif c == 1:
                walker['x'] -= 1
            elif c == 2:
                walker['y'] += 1
            else:
                walker['y'] -= 1

            x, y = walker['x'], walker['y']
            walker['d'] = min([abs(x - seed['x']) + abs(y - seed['y']) for seed in self.seeds])

        for n,walker in enumerate(self.walkers):
            if walker['d'] == 1:
                self.seeds.append({'x':walker['x'], 'y':walker['y']})
                self.walkers.pop(n)

N = 8192
L = 50
a = DLA(N, L)
while len(a.walkers) > N/2:
    a.evolve()
    print("\r\r", len(a.walkers), end="")

fig = plt.figure(figsize=(8,8))
plt.plot([p['x'] for p in a.seeds], [p['y'] for p in a.seeds], '.')
plt.show()
