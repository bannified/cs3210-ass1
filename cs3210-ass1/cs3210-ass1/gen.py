from random import *

PARTICLE_COUNTS = [100, 200, 400, 800, 1600, 3200, 6400, 12800]
L = 20000
r = 1
steps = 100
mode = "perf"

for n in PARTICLE_COUNTS:
    with open('input-%d' % n, 'w') as f:
        print(n, file=f)
        print(L, file=f)
        print(r, file=f)
        print(steps, file=f)
        print(mode, file=f)
        for i in range(n):
            print(i,
                  uniform(r, L-r),
                  uniform(r, L-r),
                  (-randint(0,1)|1) * uniform(L/4, L/8/r),
                  (-randint(0,1)|1) * uniform(L/4, L/8/r),
                  file=f)
