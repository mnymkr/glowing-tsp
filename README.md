# Travelling Salesman Problem
This project implements the Christofides algorithm in solving the Travelling Salesman Problem. It does so through these steps:
 1. Create a minimum spanning tree of the given graph.
 2. Find the set of vertices having odd degree in the minimum spanning tree.
 3. Find a minimum-weight perfect matching in the induced subgraph given by the vertices with odd degrees.
 4. Form a connected multigraph by merging the matching with the mst, each vertex how has even degree.
 5. Finally, find an Eulerian path in the multigraph, then remove the repeated vertices.
 6. We have a Hamiltonian tour.

# Requirement:
 - python 3.9

# How to run:
 1. Start by cloning the project to local.
 2. Execute christofide.py with "python christofide.py"
 3. Enter the latitude and longitude in pairs.
 4. Enter 'done' when done, the program will return a Google Maps link that contains a Hamiltonian route through every pair of coordinates.
