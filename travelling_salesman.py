from unittest import result


def hamiltonian_path(data):
    '''
    Main function for finding a Hamiltonian path in the supplied set of vertices.
    What is assumed:
        -  The graph is connected.
        -  The graph satisfies the Triangle Inequality.

    Steps to cover to produce a Hamiltonian path:
        1. Build a graph that contains all vertices from the data provided.
        2. Find a minimum spanning tree from the graph, which is a list of edges.
        4. Find all odd vertices in the minimum spanning tree, this is a list of vertices' index in the graph.
        5. Form a minimum perfect matching of all the odd vertices, which, in essence, is a list of edges.
        6. Combine the minimum spanning tree with the new minimum perfect matching to form a new
        connected multigraph, in which every vertex has even degree.
        7. Find an Eulerian tour in this new multigraph, which is a list of edges in traversal order.
        8. Finally, remove duplicate vertices from the Eulerian tour, and we'll have a Hamiltonian tour.
    '''
    # build a graph
    graph = make_graph(data)
    # print("Graph: ", G)#debug

    # build a minimum spanning tree
    minimum_spanning_tree = make_mst(graph)
    # print("Minimum Spanning Tree: ", minimum_spanning_tree)#debug

    # # find odd vertexes
    odd_vertices = find_odd_vertices(minimum_spanning_tree)
    # print("Odd vertexes in MSTree: ", odd_vertices)#debug

    # # add minimum weight matching edges to MST
    connected_multigraph = min_weight_matching(minimum_spanning_tree, graph, odd_vertices)
    # MSTree is now a connected multigraph
    # print("Connected multigraph: ", minimum_spanning_tree)

    # connected_multigraph = minimum_spanning_tree

    # find an Eulerian tour
    eulerian_tour = find_eulerian_tour(connected_multigraph, graph)
    # print("Eulerian tour: ", eulerian_tour)

    # remove duplicated nodes from the path
    current_vertex = eulerian_tour[0]
    path = [current_vertex]

    # keep a list of bools to see if one node is visited or not
    visited = [False] * len(eulerian_tour) 
    visited[eulerian_tour[0]] = True
    length = 0

    for v in eulerian_tour:

        # if v is not visited, addpend v to the path
        if not visited[v]:
            path.append(v)
            visited[v] = True

            # increment length by 1
            length += graph[current_vertex][v]

            # set current node to v
            current_vertex = v

    length += graph[current_vertex][eulerian_tour[0]]

    # first node is the last node since this is a tour
    path.append(eulerian_tour[0])

    # print("Resulting path: ", path)
    # print("Length of the resulting path: ", length)


    return path


def compute_length(x1, y1, x2, y2):
    """
    Compute the distance between 2 given nodes by their coordinates:
    x(x1,x2)
    y(y1,y2)
    
    @param 2 nodes' coordinates
    @return distance between the nodes
    """
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1.0 / 2.0)


def make_graph(data):
    """
    Continually build a graph out of the supplied data.
    The resulting graph is a dictionary
        - key: a node's index in the data list (x)
        - value: a dictionary
            + key: a node's index that (x) connects to (y)
            + value: distance between (x) and (y) computed by get_length()
    """
    graph = {}
    for this_point in range(len(data)):
        for another_point in range(len(data)):
            if this_point != another_point:
                if this_point not in graph:
                    graph[this_point] = {}

                graph[this_point][another_point] = compute_length(data[this_point][0], data[this_point][1], data[another_point][0],
                                                        data[another_point][1])

    return graph


class UnionFind:
    """
    Class of object used to store information about clusters and weight of each clusters.
    """
    def __init__(self):
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        '''
        Get the parent of the supplied object
        if the parent is not yet initialized, the object is its own parent
        otherwise, return the parent of the object'''
        # print("inside getitem") #debug
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        '''
        what
        '''
        print("iter/parents: " + self.parents)
        return iter(self.parents)

    def unionize(self, *objects):
        '''
        unionize 2 nodes into the same cluster,
        compare 2 root nodes of the clusters
        the root node with larger degree will be the other root node's parent
        '''
        roots = [self[x] for x in objects]
        x = [(self.weights[r], r) for r in roots]
        y = max(x)
        # get the heavier node's index
        heaviest = y[1]

        # the heavier root will become the other's root
        for root in roots:
            if root != heaviest:
                self.weights[heaviest] += self.weights[root]
                self.parents[root] = heaviest


def make_mst(graph):
    """
    Build a minimum spanning tree out of the supplied graph, using UnionFind 
    algorithm.

    param: 
        G - the dictionary containing the graph
    return:
        tree - the minimum spanning tree produced from the graph
    """

    tree = []
    subtrees = UnionFind()

    # sort the edges by their increasing weights
    sorted_graph = sorted((graph[u][v], u, v) for u in graph for v in graph[u])

    # go through each edge once, x and y are the vertices in that edge
    for weight, u, v in sorted_graph:

        # find x and y's root node
        root_of_u = subtrees[u]
        root_of_v = subtrees[v]

        # print("Subtree[" + str(u) + "]: " + str(sub_u))#debug
        # print("Subtree[" + str(v) + "]: " + str(sub_v))#debug

        if root_of_u != root_of_v: # x and y do not belong to the same root node

            # print("yeh not equal, unioning")#debug

            tree.append((u, v, weight)) # add the edge to the MST
            subtrees.unionize(u, v) # uninize the 2 clusters into 1

        # print()#debug
    return tree


def find_odd_vertices(minimum_spanning_tree):
    '''
    Returns a list of odd-weighted vertices in the minimum_spanning_tree.
    Build a temporary graph out of the minimum_spanning_tree for easy counting of 
    each vertex's degree.
    '''
    temp_graph = {} # used to store the degree of the vertices
    odd_vertices = []
    for edge in minimum_spanning_tree:
        # edge[0] is the first vertex incident to the edge
        if edge[0] not in temp_graph: # if not present in tmp_g, add to tmp_g
            temp_graph[edge[0]] = 0
        # edge[1] is the second vertex incident to the edge
        if edge[1] not in temp_graph:
            temp_graph[edge[1]] = 0

        # increase the degree of each vertex by 1
        temp_graph[edge[0]] += 1
        temp_graph[edge[1]] += 1

    # check if each vertex has odd degree 
    for vertex in temp_graph:
        if temp_graph[vertex] % 2 == 1: # if odd, add to odd_vetices
            odd_vertices.append(vertex)

    # print('odd vertices: ', odd_vertices)#debug
    return odd_vertices


def min_weight_matching(minimum_spanning_tree, graph, odd_vertices):
    '''
    Find a minimum weight matching in the set of odd_vertices 
    and then add the newfound vertices to the existing Minimum Spanning Tree.

    This step is supposed to be done in 2, but is instead combined since 
    new edges in the matching can be added immediately to the MST without
    any further contemplation.
    '''

    import random
    random.shuffle(odd_vertices)

    while odd_vertices: # continue until there is no odd vertices left

        # take a vertex out of odd_vertices
        current_vertex = odd_vertices.pop()

        # initial value for comparison with floating numbers
        # https://stackoverflow.com/questions/34264710/what-is-the-point-of-floatinf-in-python
        smallest_length = float("inf")

        u = 1
        closest_vertex = 0

        # begin matching every u in odd_vertices with current_vertex
        for u in odd_vertices:

            # find the shortest path between u and current_vertex
            if (current_vertex != u and # u is not current_vertex
                    graph[current_vertex][u] < smallest_length):
                smallest_length = graph[current_vertex][u]
                closest_vertex = u

        # after finding one, add the new edge to the Minimum Spanning Tree
        # duplicate edge is allowed
        minimum_spanning_tree.append((current_vertex, closest_vertex, smallest_length))

        # remove u out of odd_vertices to prevent duplication in later considerations
        odd_vertices.remove(closest_vertex)

    return minimum_spanning_tree


def find_eulerian_tour(connected_multigraph, G):
    '''
    Find an Elerian tour in the Multigraph: 
    1. Start with a vertex v
    2. Form a cycle with v.
    3. If there are any vertex (v) in that Cycle that still has edges not in that cycle yet, then:
        create a new cycle from v using unused edges,
        add that new cycle to the previous cycle, forming a larger cycle.
    4. The result will be an Eulerian cycle in the connected_graph.

    Loops of multiedges will be avoided because edges created from the perfect
    matching will only be considered after all old, 'proper' edges have been exhausted.
    '''

    # FIND NEIGHBOURS of all vertices
    neighbours = {}

    # consider each edge in MatchedMSTree
    for edge in connected_multigraph:

        # print('MatchedMSTree: ', MatchedMSTree)#debug
        
        # edge[0] is the first vertex incident to the edge
        if edge[0] not in neighbours: # if not present in neighbors, add to neighbours
            neighbours[edge[0]] = []

        # edge[1] is the second vertex incident to the edge
        if edge[1] not in neighbours:
            neighbours[edge[1]] = []

        # add each vertex in the edge as each other's neighbour
        neighbours[edge[0]].append(edge[1])
        neighbours[edge[1]].append(edge[0])
        # now neighbor looks like a dictionary: neighbors{edge[0]:[edge[1]]}
        #  - key: edge[0]
        #  - value: a list containing edge[1]

    # finds the Hamiltonian circuit
    starting_vertex = connected_multigraph[0][0] # first vertex in the first edge in connected_multigraph

    eulerian_path = [neighbours[starting_vertex][0]] # that vertex's neighbor recorded in neighbours

    # print('Current path: ', eulerian_path)#debug

    # continue until there is no edge left in connected_multigraph
    while len(connected_multigraph) > 0:

        # print('Neighbor: ', neighbours)#debug


        for index, vertex in enumerate(eulerian_path): 

            # if this vertext still has some neightbors left, then 
            # break out of the for loop with its index and vertex
            if len(neighbours[vertex]) > 0: 
                # print(str(vertex) + ' still has neighbours')#debug
                break
        
        # continue when there is still unchecked neighbour 
        # in the list of neighbours of v
        while len(neighbours[vertex]) > 0:

            # w is v's neighbour currently being considered
            neighbor = neighbours[vertex][0] 

            # remove this edge (w, v) from the connected_multigraph
            remove_edge(connected_multigraph, vertex, neighbor)

            # delete w from v's list of neighbours
            del neighbours[vertex][(neighbours[vertex].index(neighbor))]
            # delete v from w's list of neighbours
            del neighbours[neighbor][(neighbours[neighbor].index(vertex))]

            # increment index of current item in EP
            index += 1

            # insert w as the next step in the Eulerian Path
            eulerian_path.insert(index, neighbor)

            # print('Neigbor in small while: ', neighbours)#debug
            # print('Current path: ', eulerian_path)

            # start finding the next vertex, this time from w
            vertex = neighbor





    return eulerian_path


def remove_edge(connected_multigraph, v1, v2):

    '''
    Remove an edge (v1, v2) from the supplied connected_multigraph
    '''

    for i, item in enumerate(connected_multigraph):
        if (item[0] == v2 and item[1] == v1) or (item[0] == v1 and item[1] == v2):
            del connected_multigraph[i]

    return connected_multigraph

def get_map_url(result_path, mapping):

    '''
    Return a useful Google Map url which list all desired destinations in the order 
    of the Hamiltonian path of them.

    result_path is a list of indexes of items in mapping
    mapping is a list of lists, [[lat1, lon1], [lat2, lon2]]
    '''

    url = ''
    if result_path != []:
        url += 'https://www.google.com/maps/dir/?api=1'

        destination = result_path.pop()

        origin = result_path.pop(0)
        url += '&origin=' + str(mapping[origin][0]) + ',' + str(mapping[origin][1])

        url += '&mapaction=map'

        url += '&waypoints='

        for node in result_path:
            url += str(mapping[node][0]) + ',' + str(mapping[node][1]) + '|'

        url = url[:-1]

        url += '&destination=' + str(mapping[destination][0]) + ',' + str(mapping[destination][1])
    
    return url


def get_input(coordinates):

    user_input = ''
    
    while user_input != 'done':

        user_input = input("Enter lattitude and longitude, formatted as: lat,long. \nEnter 'done' when finished: ")

        if user_input == 'done':
            break

        raw_list = user_input.split(',')
        coordinates.append([float(raw_list[0]), float(raw_list[1])])

def main():
    coordinates = []
    get_input(coordinates)

    if coordinates == []:
        print("Here's your link: https://www.youtube.com/watch?v=dQw4w9WgXcQ")
        exit()

    result_path = hamiltonian_path(coordinates)

    url = get_map_url(result_path, coordinates)

    print("Here's your Google Maps link: " , url)

    # destinations = [
    #     [10.725906204249258, 106.72433446825471], # Waterfront
    #     [10.7270230306425, 106.72027991899324], #Fulbright
    #     [10.735723087013758, 106.72688522969952] # Docklands
    # ]

main()

# demo_graph = [[1, 1], [2, 5], [8, 0]]