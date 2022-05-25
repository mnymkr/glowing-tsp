from unittest import result


def hamiltonian_path(data):
    # build a graph
    graph = make_graph(data)
    # print("Graph: ", G)#debug

    # build a minimum spanning tree
    minimum_spanning_tree = minimum_spanning_tree(graph)
    # print("Minimum Spanning Tree: ", MSTree)#debug

    # # find odd vertexes
    odd_vertices = find_odd_vertices(minimum_spanning_tree)
    # print("Odd vertexes in MSTree: ", odd_vertices)

    # # add minimum weight matching edges to MST
    min_weight_matching(minimum_spanning_tree, graph, odd_vertices)
    # MSTree is now a connected multigraph
    # print("Connected multigraph: ", MSTree)

    # find an Eulerian tour
    eulerian_tour = find_eulerian_tour(minimum_spanning_tree, graph)
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

    print("Resulting path: ", path)
    print("Length of the resulting path: ", length)

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
    continually build a graph out of the supplied data.
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
    Class of object used to store information about clusters and weight of each clusters"""
    def __init__(self):
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        '''
        get the parent of the supplied object
        
        if the parent is not yet initialized, the object is its own parent
        
        else, return the parent of the object'''
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
        G the dictionary containing the graph
    return:
        tree the minimum spanning tree produced from the graph
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
    '''
    temp_graph = {}
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
    Find a minumum weight matching in the set of odd_vertices and 
    add the newfound vertices to the existing Minimum Spanning Tree.
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


def find_eulerian_tour(multigraph, G):

    # FIND NEIGHBOURS
    neighbours = {}

    # consider each edge in MatchedMSTree
    for edge in multigraph:

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
    starting_vertex = multigraph[0][0] # first vertex in the first edge in MatchedMSTree
    # EP = Eulerian Path
    eulerian_path = [neighbours[starting_vertex][0]] # that vertex's neighbor recorded in neighbours

    # print('EP: ', EP)#debug

    # continue until there is no edge left in MatchedMSTree
    while len(multigraph) > 0:
        # print('EP: ', EP)#debug
        for index, vertex in enumerate(eulerian_path): # 
            # x = neighbours[v]#debug
            if len(neighbours[vertex]) > 0:
                break # out of current for loop
        
        # continue when there is still unchecked neighbour 
        # in the list of neighbours of v
        while len(neighbours[vertex]) > 0:

            # w is v's neighbour currently being considered
            neighbor = neighbours[vertex][0] 

            # remove this edge (w, v) from the MatchedMST
            remove_edge(multigraph, vertex, neighbor)

            # delete w from v's list of neighbours
            del neighbours[vertex][(neighbours[vertex].index(neighbor))]
            # delete v from w's list of neighbours
            del neighbours[neighbor][(neighbours[neighbor].index(vertex))]

            # increment index of current item in EP
            index += 1

            # insert w as the next step in EP
            eulerian_path.insert(index, neighbor)

            # start finding the next vertex, this time from w
            vertex = neighbor

    return eulerian_path


def remove_edge(MatchedMST, v1, v2):

    for i, item in enumerate(MatchedMST):
        if (item[0] == v2 and item[1] == v1) or (item[0] == v1 and item[1] == v2):
            del MatchedMST[i]

    return MatchedMST

def get_map_url(result_path, mapping):

    # result_path is a list of indexes of items in mapping
    # mapping is a list of lists, [[lat1, lon1], [lat2, lon2]]
    url = ''
    if result_path != []:
        url += 'https://www.google.com/maps/dir/?api=1'

        # result_path.pop()
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

# each entry in tsp is [x, y]
# no need to store the length of edges because the graph is connected. Length will be computed later.
# after having the Eulerian tour, map each node's index to the place's name to make it cooler.


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
        print('No stop given, exiting...')
        exit()

    # print('coordinates: ', coordinates)

    result_path = hamiltonian_path(coordinates)

    url = get_map_url(result_path, coordinates)

    print("Here's your Google Maps link: " , url)

    # destinations = [
    #     [10.725906204249258, 106.72433446825471], # Waterfront
    #     [10.7270230306425, 106.72027991899324], #Fulbright
    #     [10.735723087013758, 106.72688522969952] # Docklands
    # ]


main()
        
# tsp([
#     [0, 0],
#     [3, 0],
#     [6, 0],

#     [0, 3],
#     [3, 3],
#     [6, 3],

#     [0, 6],
#     [3, 6],
#     [6, 6],

# ])

# tsp([[1, 1], [2, 5], [8, 0]])