# Travelling Salesman Problem
This project implements the Christofides algorithm in solving the Travelling Salesman Problem. It does so through these steps:
 1. Build a graph that contains all vertices from the data provided.
 2. Find a minimum spanning tree from the graph, which is a list of edges.
 4. Find all odd vertices in the minimum spanning tree, this is a list of vertices' index in the graph.
 5. Form a minimum perfect matching of all the odd vertices, which, in essence, is a list of edges.
 6. Combine the minimum spanning tree with the new minimum perfect matching to form a new connected multigraph, in which every vertex has even degree.
 7. Find an Eulerian tour in this new multigraph, which is a list of edges in traversal order.
 8. Finally, remove duplicate vertices from the Eulerian tour, and we'll have a Hamiltonian tour.

# Requirement:
 - python 3.9

# How to run:
 1. Start by cloning the project to local.
 2. Execute travelling_salesman.py with "python travelling_salesman.py". Alternatively, open and execute the file in your IDE.
