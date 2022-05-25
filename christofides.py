def tsp(data):
    # build a graph
    G = build_graph(data)
    print("Graph: ", G)

    # build a minimum spanning tree
    MSTree = minimum_spanning_tree(G)
    print("MSTree: ", MSTree)

    # # find odd vertexes
    odd_vertices = find_odd_vertices(MSTree)
    print("Odd vertexes in MSTree: ", odd_vertices)

    # # add minimum weight matching edges to MST
    minimum_weight_matching(MSTree, G, odd_vertices)
    print("Minimum weight matching: ", MSTree)

    # MSTree is now MatchedMSTree
    # # find an eulerian tour
    eulerian_tour = find_eulerian_tour(MSTree, G)
    print("Eulerian tour: ", eulerian_tour)

    # current = eulerian_tour[0]
    # path = [current]
    # visited = [False] * len(eulerian_tour)
    # visited[eulerian_tour[0]] = True
    # length = 0

    # for v in eulerian_tour:
    #     if not visited[v]:
    #         path.append(v)
    #         visited[v] = True

    #         length += G[current][v]
    #         current = v

    # length +=G[current][eulerian_tour[0]]
    # path.append(eulerian_tour[0])

    # print("Result path: ", path)
    # print("Result length of the path: ", length)

    # return length, path


def get_length(x1, y1, x2, y2):
    """
    Compute the distance between 2 given nodes by their coordinates:
    x(x1,x2)
    y(y1,y2)
    
    @param 2 nodes' coordinates
    @return distance between the nodes
    """
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1.0 / 2.0)


def build_graph(data):
    """
    continually build a graph out of the supplied data.
    The resulting graph is a dictionary
        - key: a node's index in the data list (x)
        - value: a dictionary
            + key: a node's index that (x) connects to (y)
            + value: distance between (x) and (y) computed by get_length()
    """
    graph = {}
    for this in range(len(data)):
        for another_point in range(len(data)):
            if this != another_point:
                if this not in graph:
                    graph[this] = {}

                graph[this][another_point] = get_length(data[this][0], data[this][1], data[another_point][0],
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

    def union(self, *objects):
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
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest


def minimum_spanning_tree(G):
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
    sorted_graph = sorted((G[u][v], u, v) for u in G for v in G[u])
    # print('sorted_graph: ', sorted_graph)#debug

    # go through each edge once, x and y are the vertices in that edge
    for Weight, u, v in sorted_graph:

        # find x and y's root node
        sub_u = subtrees[u]
        sub_v = subtrees[v]

        # print("Subtree[" + str(u) + "]: " + str(sub_u))#debug
        # print("Subtree[" + str(v) + "]: " + str(sub_v))#debug

        if sub_u != sub_v: # x and y do not belong to the same root node

            # print("yeh not equal, unioning")#debug

            tree.append((u, v, Weight)) # add the edge to the MST
            subtrees.union(u, v) # uninize the 2 clusters into 1

        # print()#debug
    return tree


def find_odd_vertices(MST):
    '''
    '''
    tmp_g = {}
    odd_vertices = []
    for edge in MST:
        # edge[0] is the first vertex incident to the edge
        if edge[0] not in tmp_g: # if not present in tmp_g, add to tmp_g
            tmp_g[edge[0]] = 0
        # edge[1] is the second vertex incident to the edge
        if edge[1] not in tmp_g:
            tmp_g[edge[1]] = 0

        # increase the degree of each vertex by 1
        tmp_g[edge[0]] += 1
        tmp_g[edge[1]] += 1

    # check if each vertex has odd degree 
    for vertex in tmp_g:
        if tmp_g[vertex] % 2 == 1: # if odd, add to odd_veti
            odd_vertices.append(vertex)

    print('odd vertices: ', odd_vertices)#debug
    return odd_vertices


def minimum_weight_matching(MST, G, odd_vertices):
    '''
    Find a minumum weight matching in the set of odd_vertices and 
    add the newfound vertices to the existing Minimum Spanning Tree.
    '''

    import random
    random.shuffle(odd_vertices)

    while odd_vertices: # continue until there is no odd vertices left

        # take v out of odd_vertices
        v = odd_vertices.pop()

        # initial value for comparison with floating numbers
        # https://stackoverflow.com/questions/34264710/what-is-the-point-of-floatinf-in-python
        length = float("inf")

        u = 1
        closest = 0

        # begin matching every u in odd_vertices with v
        for u in odd_vertices:

            # find the shortest path between v and v
            if v != u and G[v][u] < length:
                length = G[v][u]
                closest = u

        # after finding one, add the new edge to the Minimum Spanning Tree
        # duplicate edges do not matter
        MST.append((v, closest, length))

        # remove u out of odd_vertice to prevent duplication in later considerations
        odd_vertices.remove(closest)


def find_eulerian_tour(MatchedMSTree, G):

    # find neigbours
    neighbours = {}

    # consider each edge in MatchedMSTree
    for edge in MatchedMSTree:
        
        # edge[0] is the first vertex incident to the edge
        if edge[0] not in neighbours: # if not present in neighbors, add to neighbours
            neighbours[edge[0]] = []

        # edge[1] is the second vertex indident to the edge
        if edge[1] not in neighbours:
            neighbours[edge[1]] = []

        # add each vertex in the edge as each other's neighbour
        neighbours[edge[0]].append(edge[1])
        neighbours[edge[1]].append(edge[0])
        # now neighbor looks like a dictionary: neighbors{edge[0]:[edge[1]]}
        #  - key: edge[0]
        #  - value: a list containing edge[1]

    # print("Neighbours: ", neighbours)

    # finds the hamiltonian circuit
    start_vertex = MatchedMSTree[0][0] # first vertex in the first edge in MatchedMSTree
    EP = [neighbours[start_vertex][0]] # 

    # continue until there is no edge left in MatchedMSTree
    while len(MatchedMSTree) > 0:
        for i, v in enumerate(EP):
            if len(neighbours[v]) > 0:
                break

        while len(neighbours[v]) > 0:
            w = neighbours[v][0]

            remove_edge_from_matchedMST(MatchedMSTree, v, w)

            del neighbours[v][(neighbours[v].index(w))]
            del neighbours[w][(neighbours[w].index(v))]

            i += 1
            EP.insert(i, w)

            v = w

    return EP


def remove_edge_from_matchedMST(MatchedMST, v1, v2):

    for i, item in enumerate(MatchedMST):
        if (item[0] == v2 and item[1] == v1) or (item[0] == v1 and item[1] == v2):
            del MatchedMST[i]

    return MatchedMST


def print_google_map():
    """This function prints out a link to Google Maps, where you can view ISS movements
    during the time which function track_iss() was called. This function only works
    if you have called track_iss() before."""
    if long_lats != []:
        print('Go here to see a Google Map of path followed by ISS: https://www.google.com/maps/dir', end='')
        for i in range(len(long_lats)):
            # This loop adds all coordinations to the link
            print('/'+long_lats[i], end='')
        print('//@'+long_lats[len(long_lats)//2]+',5z')
        #This command sets the center of the map to the mid-point of ISS trajactory,
        #it also sets the zoom level to 5z, enough to view the trajactory in detail.
    else:
        print('Please run track_iss BEFORE running this function')7

# each entry in tsp is [x, y]
# no need to store the length of edges because the graph is connected. Length will be computed later.
# after having the Eulerian tour, map each node's index to the place's name to make it cooler.

# tsp([[1380, 939], [2848, 96], [3510, 1671], [457, 334], [3888, 666], [984, 965], [2721, 1482], [1286, 525],
#              [2716, 1432], [738, 1325], [1251, 1832], [2728, 1698], [3815, 169], [3683, 1533], [1247, 1945], [123, 862],
#              [1234, 1946], [252, 1240], [611, 673], [2576, 1676], [928, 1700], [53, 857], [1807, 1711], [274, 1420],
#              [2574, 946], [178, 24], [2678, 1825], [1795, 962], [3384, 1498], [3520, 1079], [1256, 61], [1424, 1728],
#              [3913, 192], [3085, 1528], [2573, 1969], [463, 1670], [3875, 598], [298, 1513], [3479, 821], [2542, 236],
#              [3955, 1743], [1323, 280], [3447, 1830], [2936, 337], [1621, 1830], [3373, 1646], [1393, 1368],
#              [3874, 1318], [938, 955], [3022, 474], [2482, 1183], [3854, 923], [376, 825], [2519, 135], [2945, 1622],
#              [953, 268], [2628, 1479], [2097, 981], [890, 1846], [2139, 1806], [2421, 1007], [2290, 1810], [1115, 1052],
#              [2588, 302], [327, 265], [241, 341], [1917, 687], [2991, 792], [2573, 599], [19, 674], [3911, 1673],
#              [872, 1559], [2863, 558], [929, 1766], [839, 620], [3893, 102], [2178, 1619], [3822, 899], [378, 1048],
#              [1178, 100], [2599, 901], [3416, 143], [2961, 1605], [611, 1384], [3113, 885], [2597, 1830], [2586, 1286],
#              [161, 906], [1429, 134], [742, 1025], [1625, 1651], [1187, 706], [1787, 1009], [22, 987], [3640, 43],
#              [3756, 882], [776, 392], [1724, 1642], [198, 1810], [3950, 1558]])

tsp([[1, 1], [2, 5], [8, 0]])

#
# tsp([
#     [0, 0],
#     [3, 0],
#     [6, 0],
#
#     [0, 3],
#     [3, 3],
#     [6, 3],
#
#     [0, 6],
#     [3, 6],
#     [6, 6],
#
# ])
