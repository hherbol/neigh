# neigh
Extension for Python - Generate neighbour list

## Description
This code was written to allow for the faster generation of neighbour lists, given a list of elements in n-dimensions.
It uses a supplied cutoff radius, and currently brute-forces the calculation.  Further iterations will allow for more
efficient methods.

## How to use
Simply place the *.so* or *.pyd* file in the same directory as your code, depending on if you're using linux or windows
respectively, then import.

    import neigh
    from random import random

    X, Y, Z = 100, 100, 100
    O_x, O_y, O_z = 0, 0, 0
    N = 100
    CUTOFF = 1.5

    POSITION_ARRAY = [None for i in range(N)]
    for i in range(N):
        x, y, z = random()*X+O_x, random()*Y+O_y, random()*Z+O_z
        positions[i] = [x, y, z]

    neighbour_list = neigh.get_neigh(POSITION_ARRAY, CUTOFF, PBC=[X, Y, Z], origin = [O_x, O_y, O_z])

The function `neigh.get_neigh` takes in:

- POSITION_ARRAY = A list of lists, each holding floats describing your system
- CUTOFF = A float that gives a cutoff radius
- PBC = A list of floats describing your box
- origin = A list of floats describing where your origin is

NOTE! In theory the code allows for n-dimensions; however, if you use PBC, then it only allows for dim <= 3.
