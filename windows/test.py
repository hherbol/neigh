from neigh import gen_neigh
from random import random # Rand 0.0 <= val < 1.0
import time

X, Y, Z = 1000,1000,1000
o1, o2, o3 = 0, 0, 0
N = 8000
CUTOFF = 10.0

positions = [None for i in range(N)]
for i in range(N):
	x, y, z = random()*X+o1, random()*Y+o2, random()*Z+o3
	positions[i] = [x,y,z]

print("\nStarting..."),
t0 = time.clock()
neighbours = gen_neigh(positions, CUTOFF, PBC=[X,Y,Z], origin=[o1,o2,o3])
print(" Done\n")

print("Time taken = %.6f" % (time.clock()-t0))

