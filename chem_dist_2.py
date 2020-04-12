#!/usr/bin/python3

from collections import defaultdict
import copy
import operator
import numpy as np
import pygame

inf = 10**10

### Operations on the plane ###

def tuplesum(x,y):
    return tuple(map(operator.add,x,y))

def tuplescalarprod(c,x):
    return tuple([c*coord for coord in x])

# def inbox(x,n):
#     for coord in x:
#         if coord < 0 or coord > n:
#             return False
#     return True

def nbrs(x,n):
    xnbrs = []
    direction = [(-2,0),(2,0),(0,-2),(0,2)]
    for y in direction:
        nbr = tuplesum(x,y)
        if 0 <= nbr[0] and nbr[0] <= n and 0 <= nbr[1] and nbr[1] <= n:
            xnbrs.append(nbr)
        # if inbox(alt,n):
        #     xnbr.append(alt)
    return xnbrs

### edge append for undirected graph

def edge_append(x,y,w,G):
    G[x][y] = w
    G[y][x] = w
    return G

### Generate a configuration of a percolation with parameter p on a box of size n

def ber(p):
    if np.random.random() <= p:
        return 1
    else:
        return inf

def perc(n,p):
    G = defaultdict(dict)
    for x1 in range(0,n+2,2):
        for x2 in range(0,n+2,2):
            #x = tuplescalarprod(2,(x1,x2))
            x = (x1,x2)
            xnbrs = nbrs(x,n)
            for y in xnbrs:
                if y not in list(G[x].keys()):
                    G = edge_append(x,y,ber(p),G)
    return G

### Dual graph

def edge_dual(x,y):
    e = [x,y]
    e.sort()
    if e[0][0] == e[1][0]: # e is vertical
        e = [tuplesum(e[0],(1,1)),tuplesum(e[1],(-1,-1))]
    elif e[0][1] == e[1][1]: # e is horizontal
        e = [tuplesum(e[0],(1,-1)),tuplesum(e[1],(-1,1))]
    return e

def dual_weight(w):
    if w == 1:
        return inf
    elif w == inf:
        return 1

def dual(G,n):
    G_dual = defaultdict(dict)
    vert = G.keys()
    for x in vert:
        xnbrs = nbrs(x,n)
        for y in xnbrs:
            [x_dual,y_dual] = edge_dual(x,y)
            edge_append(x_dual,y_dual,dual_weight(G[x][y]),G_dual)
    return G_dual

### Graph drawer

def graph_draw(screen,G,n,edge_length,color,width):
    #fig_size = edge_length*n
    #screen = pygame.Surface((fig_size,fig_size))
    vert = G.keys()
    for x in vert:
        xnbrs = list(G[x].keys())
        for y in xnbrs:
            if G[x][y] == 1:
                xcoord = tuplescalarprod(edge_length,x)
                ycoord = tuplescalarprod(edge_length,y)
                pygame.draw.line(screen,color,xcoord,ycoord,width)
    return screen

def path_draw(screen,path,edge_length,color,width):
    path = [tuplescalarprod(edge_length,x) for x in path]
    pygame.draw.lines(screen,color,False,path,width)

### Path existence

def connected(G, start, end, visited=None):
    if start == end:
        return True
    if visited is None:
        visited = []
    visited.append(start)
    found = False
    nbr = list(G[start].keys())
    nbr_not_visited = []
    for item in nbr:
        if item not in visited and G[start][item] < inf:
            nbr_not_visited.append(item)
    nbr_not_visited = [item for item in nbr if item not in visited and G[start][item] < inf]
    for node in nbr_not_visited:
        if found:
            return True
        if node == end:
            return True
        found = connected(G, node, end, visited)
    return found

### Shortest path

def spath(G,start,end):
    [mindist,prevert] = dijkstra(G,start)
    spath_list = [end]
    while prevert[end] != start:
        end = prevert[end]
        spath_list.append(end)
    spath_list.append(start)
    return spath_list

# python will change the input variable
def dijkstra(G,start):
    unvisited = list(G.keys())
    mindist = {}
    for x in unvisited:
        mindist[x] = inf
    mindist[start] = 0
    prevert = {}
    while unvisited:
        v = unvisited[0]
        for u in unvisited:
            if mindist[u] < mindist[v]:
                v = u
        vnbr = list(G[v].keys())
        for w in vnbr:
            alt = mindist[v]+G[v][w]
            if mindist[w] > alt:
                mindist[w] = alt
                prevert[w] = v
        unvisited.remove(v)
    return [mindist,prevert]

### Cluster exploration: return the cluster on G containing x
def cluster(x,G):
    C = defaultdict(dict)
    C_bdy = defaultdict(dict)
    C[x] = G[x]
    visited = [x]
    queue = [x]
    while queue:
        y = queue.pop(0)
        C[y] = G[y]
        ynbrs = [z for z in list(G[y].keys()) if G[y][z] < inf]
        for nbr in ynbrs:
            if nbr not in visited:
                visited.append(nbr)
                queue.append(nbr)
    return C

### Find the edges around the cluster

def nbrs_dual(x,G_dual):
    xnbrs_dual = []
    direction_dual = [(1,1),(-1,1),(-1,-1),(1,-1)]
    for y in direction_dual:
        z = tuplesum(x,y)
        if z in G_dual:
            xnbrs_dual.append(z)
    return xnbrs_dual

def lpath(G,G_dual,bottom,start,end):
    C_dual_bdy = defaultdict(dict)
    C_dual_bdy[start] = G[start]
    C_dual_bdy[end] = G[end]
    C = cluster(bottom,G_dual)
    del C[bottom] # bottom is not on the graph
    vertC = list(C.keys())
    for x_dual in vertC:
        xnbrs_dual = nbrs_dual(x_dual,G)
        for y in xnbrs_dual:
            vertC_dual_bdy = list(C_dual_bdy.keys())
            if y not in vertC_dual_bdy:
                C_dual_bdy[y] = G[y]
    found = False
    path = []
    [found, path] = dfs(C_dual_bdy,start,end,path,found)
    return [path,C]

def dfs(G,start,end,path,found,visited = None):
    if start == end:
        path.append(start)
        return [True, path]
    if visited == None:
        visited = []
    visited.append(start)
    path.append(start)
    nbr = list(G[start].keys())
    nbr_not_visited = [item for item in nbr if item not in visited and G[start][item] < inf]
    for node in nbr_not_visited:
        if found:
            return [True, path]
        if node == end:
            path.append(node)
            return [True, path]
        [found, path] = dfs(G, node, end, path, found, visited)
    if not found:
        path = path[:-1] # Doesn't end at the main branch
    return [found, path]

# G = {
# 0: {1: 1, 4: 1},
# 1: {0: 1},
# 2: {3: 1, 5: 1},
# 3: {2: 1},
# 4: {0: 1, 5: 1},
# 5: {2: 1, 4: 1, 8: 1},
# 6: {10: 1},
# 7: {8: 1},
# 8: {5: 1, 7: 1, 9: 1},
# 9: {8: 1, 10: 1},
# 10: {6: 1, 9: 1},
# }
# start = 0
# end = 10
# print(dfs(G,start,end,[],False))

# Test example for cluster
# G = {
# 1: {2: inf},
# 2: {3: 1},
# 3: {2: 1},
# }
#
# x = 3
# C = cluster(x,G)
# if C:
#     print(C)
# else:
#     print("cluster does not exist")

### Main Part ###
### Working on 2Zx2Z grid


#print("Side length of the box")
#n = 2*int(input())
def simulation(n):
    f = open("paths_info", "a")
    print("Side length of the box = {}".format(n), file = f)
    p = 0.5
    G = perc(n,p)
    G_dual = dual(G,n)

    lG = [x for x in G.keys() if x[0] == 0]
    rG = [x for x in G.keys() if x[0] == n]
    bG_dual = [x for x in G_dual.keys() if x[1] == -1]
    uG_dual = [x for x in G_dual.keys() if x[1] == n+1]

    start = 'start'
    end = 'end'

    def start_end_add(G,START,END):
        for x in START:
            edge_append(start,x,0,G)
        for x in END:
            edge_append(end,x,0,G)
        return G

    G = start_end_add(G,lG,rG)
    G_dual = start_end_add(G_dual,bG_dual,uG_dual)

    # for x in lG:
    #     edge_append(start,x,0,G) # dist(start,lG) = 0 on G
    # for x in rG:
    #     edge_append(end,x,0,G) # dist(end,rG) = 0 on G
    # for x in bG_dual:
    #     edge_append(start,x,0,G_dual) # dist(start,bG_dual) = 0 on G_dual
    # for x in uG_dual:
    #     edge_append(end,x,0,G_dual) # dist(end,uG_dual) = 0 on G_dual


    ### The part that draws the percolation ###

    edge_length = 10
    color_prime = (220,220,220) # Gainsboro Grey
    color_dual = (169,169,169) # Dark Grey
    color_spath = (0,0,0) # Black
    color_lpath = (255,0,0) # Dark Red
    color_cluster = (220,220,220) # Dark Grey
    edge_width = 2
    path_width = 3
    cluster_width = 4
    #fig_size = edge_length*(n+12)
    fig_size = edge_length*(n)
    screen = pygame.Surface((fig_size,fig_size))
    screen.fill(pygame.Color("white"))
    # graph_draw(screen,G,n,edge_length,color_prime,edge_width)
    # pygame.image.save(screen,'./chem_dist_sym.png')
    # screen2 = graph_draw(screen,G_dual,n,edge_length,color_dual,edge_width)
    # pygame.image.save(screen2,'./chem_dist_sym_2.png')
    if connected(G,start,end):
        print("Horizontal crossing exists.")
        graph_draw(screen,G,n,edge_length,color_prime,edge_width)
        pygame.image.save(screen,'./chem_dist_sym_1.png')
        # Lowest path
        [lpath_list,C] = lpath(G,G_dual,start,start,end)
        lpath_list.remove(start)
        lpath_list.remove(end)
        print("Length of the lowest path = {}.".format(len(lpath_list)-1), file = f)
        print(lpath_list, file = f)
        path_draw(screen,lpath_list,edge_length,color_lpath,path_width)
        pygame.image.save(screen,'./chem_dist_sym_2.png')
        graph_draw(screen,C,n,edge_length,color_dual,cluster_width)
        pygame.image.save(screen,'./chem_dist_sym_3.png')
        # Shortest path
        spath_list = spath(G,start,end)
        spath_list.remove(start)
        spath_list.remove(end)
        print("Length of the shortest path = {}.".format(len(spath_list)-1), file = f)
        print(spath_list, file = f)
        path_draw(screen,spath_list,edge_length,color_spath,path_width)
        pygame.image.save(screen,'./chem_dist_sym_4.png')
    else:
        print("No horizontal crossing. Consider the dual graph.")
        graph_draw(screen,G_dual,n,edge_length,color_prime,edge_width)
        pygame.image.save(screen,'./chem_dist_sym_1.png')
        # Lowest path
        [lpath_list,C] = lpath(G_dual,G,start,start,end)
        lpath_list.remove(start)
        lpath_list.remove(end)
        print("Length of the lowest path = {}.".format(len(lpath_list)-1), file = f)
        print(lpath_list, file = f)
        path_draw(screen,lpath_list,edge_length,color_lpath,path_width)
        pygame.image.save(screen,'./chem_dist_sym_2.png')
        graph_draw(screen,C,n,edge_length,color_dual,cluster_width)
        pygame.image.save(screen,'./chem_dist_sym_3.png')
        # Shortest path
        spath_list = spath(G_dual,start,end)
        spath_list.remove(start)
        spath_list.remove(end)
        print("Length of the shortest path = {}.".format(len(spath_list)-1), file = f)
        print(spath_list, file = f)
        path_draw(screen,spath_list,edge_length,color_spath,path_width)
        pygame.image.save(screen,'./chem_dist_sym_4.png')
    f.close
    return
