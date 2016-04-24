from neigh import gen_neigh
from random import random # Rand 0.0 <= val < 1.0

LENGTH = 10.0
N = 10
CUTOFF = 5.2



positions = [None for i in range(N)]
for i in range(N):
	x, y, z = random()*LENGTH, random()*LENGTH, random()*LENGTH
	positions[i] = [x,y,z]

print("\nPOSITIONS: %s" % str(positions))

print("\n\nA..."),
neighbours = gen_neigh(positions, CUTOFF)
print("B...\n")

print("\n\nNeighs: %s" % str(neighbours))

#print("Returned %d" % neighbours)

#neighbours2 = gen_neigh(positions, R = 1, PBC=[100.0, 100.0, 100.0])

