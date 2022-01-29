# -*- coding: utf-8 -*-
"""
functions for DeCAF

Created on Wed Apr  8 09:41:59 2015
by Marta Stepniewska
"""

from decaf import PHARS, COLORS, Pharmacophore
import numpy as np
import math
import warnings

warnings.simplefilter("always", UserWarning)


def compare_nodes(n1, n2):
    """Compare types of two nodes. Return unnormalised similarity score and new
    dictionary of pharmacophoric properties for nodes combination.

    Args:
       n1, n2 (dict): nodes to compare

    Returns:
      float: unnormalised similarity score
      dict: pharmacophoric properties for nodes combination
    """
    if not isinstance(n1, dict):
        raise TypeError("Invalid n1! Expected dict, got %s instead" %
                        type(n1).__name__)
    if not isinstance(n2, dict):
        raise TypeError("Invalid n2! Expected dict, got %s instead" %
                        type(n2).__name__)

    if not Pharmacophore.check_node(n1):
        raise ValueError("Invalid n1!")

    if not Pharmacophore.check_node(n2):
        raise ValueError("Invalid n2!")

    c = n1["freq"] + n2["freq"]
    d1 = sum(n1["type"].values())
    d2 = sum(n2["type"].values())
    d = d1 + d2
    sim = 0.0
    t = {}

    for phar in PHARS:
        if phar in n1["type"] and phar in n2["type"]:
            sim += (n1["type"][phar] + n2["type"][phar]) / d
            t[phar] = n1["type"][phar] + n2["type"][phar]
        elif phar in n1["type"]:
            t[phar] = n1["type"][phar]
        elif phar in n2["type"]:
            t[phar] = n2["type"][phar]
    return sim * c, t


def get_rings(phar):
    """Find ring systems in given Pharmacophore. Note that only nodes of type
    "R" are considered.

    Args:
       phar (Pharmacophore): pharmacophore

    Returns:
       list of dicts: nodes representing compressed ring systems
       list of lists: members of found ring systems
    """

    if not isinstance(phar, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(phar).__name__)

    def dfs_backedge(p, n, to_check=None, visited=None, spanning_tree=None):

        cycles = []
        if visited is None:
            visited = []

        if to_check is None:
            to_check = set(range(p.numnodes))

        if spanning_tree is None:
            spanning_tree = {n: None}

        tmp = list(to_check)

        for v in tmp:
            if v in np.where(p.edges[n] > 0.0)[0]:
                if v not in visited:
                    visited.append(v)
                    to_check.remove(v)
                    spanning_tree[v] = n
                    cycles += dfs_backedge(p, v, to_check, visited,
                                           spanning_tree)
                elif spanning_tree[n] != v:
                    w = n
                    cycle = set([v])
                    add = True
                    while w != v:
                        v = spanning_tree[v]
                        cycle.add(v)
                    if add:
                        cycles.append(cycle)
        return cycles

    rings_members = set()
    for n in range(phar.numnodes):
        if "R" in phar.nodes[n]["type"]:
            rings_members.add(n)

    cycles = []
    while len(rings_members) > 0:
        node = rings_members.pop()
        cycles += dfs_backedge(phar, node, to_check=rings_members)

    # join fused ring systems
    to_del = []
    for i in range(len(cycles)):
        for j in range(i):
            if len(cycles[i] & cycles[j]) > 0:
                cycles[i] = (cycles[i] | cycles[j])
                to_del.append(j)

    for i in range(len(cycles)-1, -1, -1):
        if i in to_del:
            del cycles[i]
        else:
            cycles[i] = list(cycles[i])

    ring_nodes = []
    for i in range(len(cycles)):
        n = float(len(cycles[i]))
        ring_node = {"label": "R"+str(i), "freq": 0.0, "type": {}}

        for j in cycles[i]:
            ring_node["freq"] += phar.nodes[j]["freq"]
            for t in phar.nodes[j]["type"]:
                if t not in ring_node["type"]:
                    ring_node["type"][t] = phar.nodes[j]["type"][t]
                else:
                    ring_node["type"][t] += phar.nodes[j]["type"][t]

        ring_nodes.append(ring_node)

    return ring_nodes, cycles


def distances(p):
    """Compute lengths of shortest paths between all nodes in Pharmacophore.

    Args:
       p (Pharmacophore): model to analyse

    Returns:
       dist (numpy array): array with distances between all nodes
    """

    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    dist = np.array(p.edges)

    for i in range(p.numnodes):
        for j in range(i):
            if dist[i][j] == 0:
                dist[i][j] = dist[j][i] = float("inf")

    for i in range(len(dist)):
        compute = False
        for j in range(i):
            if dist[i][j] == float("inf"):
                compute = True
                break
        if compute:
            queue = [k for k in range(p.numnodes)]
            while queue:
                queue.sort(key=lambda x: dist[i, x])
                u = queue[0]
                del queue[0]
                for v in np.where(p.edges[u] > 0)[0]:
                    if v in queue:
                        alt = dist[i, u] + p.edges[u, v]
                        if alt < dist[i, v]:
                            dist[i, v] = dist[v, i] = alt
    return dist


def dfs(p, n, to_check=None, visited=None):
    """Perform depth-first search.

    Args:
       p (Pharmacophore): model to search
       n (int): id of first node
       to_check (set, optional): indices of nodes to check
       visited (list, optional): list of indicies of already visited nodes; if
       given it will be updated

    Returns:
       visited (list): all nodes reachable from n
    """
    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    if visited is None:
        visited = []
    elif not isinstance(visited, list):
        raise TypeError("Expected list, got %s instead" %
                        type(visited).__name__)
    else:
        for i in visited:
            if not isinstance(i, int):
                raise TypeError("Invalid visited list! Node %s is not int!" % i)
            elif i < 0 or i >= p.numnodes:
                raise ValueError("Invalid visited list! Node index out of "
                                 "range: %s" % i)

    if to_check is None:
        to_check = set(range(p.numnodes))
    elif not isinstance(to_check, set):
        raise TypeError("Expected set, got %s instead" % type(to_check).__name__)
    else:
        for i in to_check:
            if not isinstance(i, int):
                raise TypeError("Invalid to_check list! Node %s is not int!" % i)
            elif i < 0 or i >= p.numnodes:
                raise ValueError("Invalid to_check list! Node index out of"
                                 "range: %s" % i)
    if not isinstance(n, int):
        raise TypeError("Starting node is not int!")
    elif n < 0 or n >= p.numnodes:
        raise ValueError("Starting node index out of range!")

    tmp = list(to_check)
    for v in tmp:
        if v not in visited:
            if p.edges[n, v] > 0.0:
                visited.append(v)
                to_check.remove(v)
                dfs(p, v, to_check, visited)
    return visited


def split_components(p, nodes=None):
    """Find all connected components in given Pharmacophore.

    Args:
       p (Pharmacophore): model to analyse
       nodes (list, optional): list of nodes indices; if given, find
         components in subgraph induced by those nodes.

    Returns:
       list: nodes indices grouped into connected components, sorted by
         component size
    """
    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    if nodes is None:
        nodes = list(range(p.numnodes))
    else:
        if not isinstance(nodes, list):
            raise TypeError("Expected nodes list, got %s instead" %
                            type(nodes).__name__)
        else:
            for n in nodes:
                if not isinstance(n, int):
                    raise TypeError("Node %s is not int!" % n)
                elif n < 0 or n >= p.numnodes:
                    raise ValueError("Node index out of range: %s" % n)

    to_check = set(nodes)
    visited = []
    components = []

    for n in nodes:
        if n in to_check:
            visited = [n]
            to_check.remove(n)
            dfs(p, n, to_check, visited)
            components.append(visited[:])

    return sorted(components, key=len, reverse=True)


def __modular_product(p1, p2, dist1=None, dist2=None, dist_tol=0):
    """Create modular product of given pharmacphores, that can be used to find
    their best coarse-grained alignment. All ring systems are compressed to
    single nodes and treated as sets of features. Distances between nodes are
    used as constraints to reduce number of edges in the graph.

    Args:
       p1, p2 (Pharmacophore): models to align
       dist1, dist2 (numpy array, optional): arrays with distances between all
         nodes in models
       dist_tol (float, optional): accept distance differences below this
         threshold

    Returns:
       modular product graph, consists of:
         list of dicts: nodes
         1D numpy array: partial scores coresponging to nodes
         2D numpy array: edges
         2D numpy array: length differences costs for all edges
    """

    if dist1 is None:
        dist1 = distances(p1)
        dist1[p1.edges > 0] = p1.edges[p1.edges > 0]

    if dist2 is None:
        dist2 = distances(p2)
        dist2[p2.edges > 0] = p2.edges[p2.edges > 0]

    nodes = []
    scores = []
    rings1, members1 = get_rings(p1)
    rings2, members2 = get_rings(p2)

    rings_members1 = [node for cycle in members1 for node in cycle]
    rings_members2 = [node for cycle in members2 for node in cycle]

    for i in range(p1.numnodes):
        for j in range(p2.numnodes):
            if i in rings_members1 and j in rings_members2:
                # do not align parts of rings
                # whole rings will be aligned later
                continue
            weighted_freq, _ = compare_nodes(p1.nodes[i], p2.nodes[j])
            if weighted_freq > 0.0:
                nodes.append({"n1": i, "n2": j})
                scores.append(weighted_freq)

    for i in range(len(rings1)):
        for j in range(len(rings2)):
            # do not look at the number of nodes in the ring
            # it would be faster, but it results in wrong alignemtn for complex
            # ring systems created from multiple molecules
            weighted_freq, _ = compare_nodes(rings1[i], rings2[j])
            if weighted_freq > 0.0:
                nodes.append({"n1": p1.numnodes+i, "n2": p2.numnodes+j,
                              "members": [members1[i], members2[j]]})
                scores.append(weighted_freq)

    n = len(nodes)
    scores = np.array(scores)
    edges = np.zeros((n, n))
    costs = np.zeros((n, n))

    for i in range(n):
        for j in range(i):

            if nodes[i]["n1"] == nodes[j]["n1"] or \
               nodes[i]["n2"] == nodes[j]["n2"]:
                continue

            if nodes[i]["n1"] >= p1.numnodes:   # ring node
                u = nodes[i]["n1"] - p1.numnodes
                v = nodes[i]["n2"] - p2.numnodes
                # get all nodes forming a ring system
                idxi1 = members1[u]
                idxi2 = members2[v]

            else:
                u = nodes[i]["n1"]
                v = nodes[i]["n2"]
                idxi1 = [u]
                idxi2 = [v]

            if nodes[j]["n1"] >= p1.numnodes:   # ring node
                w = nodes[j]["n1"] - p1.numnodes
                s = nodes[j]["n2"] - p2.numnodes
                idxj1 = members1[w]
                idxj2 = members2[s]

            else:
                w = nodes[j]["n1"]
                s = nodes[j]["n2"]
                idxj1 = [w]
                idxj2 = [s]

            if len(set(idxi1) & set(idxj1)) > 0 or \
               len(set(idxi2) & set(idxj2)) > 0:
                # do not connect node with itself
                continue

            is_connected = False
            # compute distances in graphs
            # for ring nodes find shortest distances
            # note: loop is faster than numpy (a lot of singletons to check)
            d1 = float("inf")
            for p in idxi1:
                for q in idxj1:
                    if p1.edges[p, q] > 0:
                        is_connected = True
                    if dist1[p, q] < d1:
                        d1 = dist1[p, q]

            d2 = float("inf")
            for p in idxi2:
                for q in idxj2:
                    if p2.edges[p, q] > 0:
                        is_connected = True
                    if dist2[p, q] < d2:
                        d2 = dist2[p, q]

            if math.fabs(d1 - d2) <= dist_tol:
                if is_connected:
                    costs[i, j] = costs[j, i] = math.fabs(d1 - d2)
                edges[i, j] = edges[j, i] = 1.0
    return nodes, scores, edges, costs


def __BronKerbosch(edges, P=None, X=None, R=None, degrees=None, neigh=None):
    """Bron-Kerbosch algorithm for finding all maximal cliques in a graph

    Args:
       edges (numpy array): array representing edges in the graph
       P (set of ints, optional): nodes to check
       X (set of ints, optional): excluded nodes
       R (set of ints, optional): current clique
       degrees (numpy array, optional): nodes degrees
       neigh (dict, optional): dictionary of sets of neighbours for all nodes

    Returns:
       list of sets: list of all maximal cliques

    see:
       Bron C, Kerbosch J. "Algorithm 457: finding all cliques of an undirected
       graph." Commun ACM. 1973;16(9):575–577.

       Cazals F, Karande C. "A note on the problem of reporting maximal
       cliques." Theor Comput Sci. 2008;407(1–3):564–568.
    """

    if P is None:
        P = set(range(len(edges)))
    if X is None:
        X = set()
    if R is None:
        R = set()

    if degrees is None:
        degrees = np.sum(edges > 0, axis=1)

    if neigh is None:
        neigh = {}
        for i in range(len(edges)):
            neigh[i] = set(np.where(edges[i] > 0)[0])

    if len(P) == 0 and len(X) == 0:
        yield R
    else:
        candidates = np.array(list(P | X))
        # try to select pivot which minimizes number of recursive calls
        pivot = candidates[np.argmax(degrees[candidates])]

        for v in (P - neigh[pivot]):
            for clique in __BronKerbosch(edges, degrees=degrees,
                                         R=(R | set([v])),
                                         P=(P & neigh[v]),
                                         X=(X & neigh[v]),
                                         neigh=neigh):
                yield clique
            P = P - set([v])
            X = X | set([v])


def __align_rings(p1, p2, n1, n2, idx1, idx2, mapping=None, dist1=None,
                  dist2=None, dist_tol=0):
    """Align rings from coarse-grained alignment.

    Args:
       p1, p2 (Pharmacophore): models to align
       n1, n2 (list of numpy arrays): nodes to align. i-th array of both
         lists contains parts of coarse-grained alignment (i.e. ring systems)
         that should be mapped to each other.
       idx1, idx2 (list of ints): lists of aligned nodes
       mapping (numpy array, optional): array describing nodes compatibility
       dist1, dist2 (numpy array, optional): arrays with distances between all
         nodes in models
       dist_tol (float, optional): accept distance differences below this
         threshold

    Returns:
       float: unnormalized similarity score
       float: edge length differences cost
       2D list: list of two lists representing matched nodes
    """

    assert len(n1) == len(n2), "wrong n1 or n2"

    if mapping is None:
        mapping = np.zeros((p1.numnodes, p2.numnodes))

        for i in range(p1.numnodes):
            for j in range(p2.numnodes):
                weighted_freq, _ = compare_nodes(p1.nodes[i], p2.nodes[j])
                if weighted_freq > 0.0:
                    mapping[i][j] = weighted_freq

    if dist1 is None:
        dist1 = distances(p1)
        dist1[p1.edges > 0] = p1.edges[p1.edges > 0]

    if dist2 is None:
        dist2 = distances(p2)
        dist2[p2.edges > 0] = p2.edges[p2.edges > 0]

    # create modular product of graph using alignment as constraints
    nodes = []
    scores = []

    assert len(idx1) == len(idx2), "unequal subgraphs sizes"
    old_len = len(idx1)

    for i in range(old_len):
        weighted_freq = mapping[idx1[i], idx2[i]]
        assert weighted_freq > 0, "wrong alignment given"
        nodes.append({"n1": idx1[i], "n2": idx2[i]})
        scores.append(weighted_freq)

    for i in range(len(n1)):
        possible_matches = np.where(mapping[n1[i], :][:, n2[i]] > 0)
        pairs1 = n1[i][possible_matches[0]]
        pairs2 = n2[i][possible_matches[1]]
        if len(idx1) > 0:
            # find pairs compatible with given alignment
            d1 = dist1[pairs1, :][:, idx1]
            d2 = dist2[pairs2, :][:, idx2]
            compatible = np.where(np.all(np.abs(d1 - d2) <= dist_tol,
                                         axis=1))[0]
        else:
            # empty alignment given, accept everything
            compatible = np.array(list(range(len(pairs1))))

        for i in compatible:
            weighted_freq = mapping[pairs1[i], pairs2[i]]
            assert weighted_freq > 0, "wrong possible matches"
            nodes.append({"n1": pairs1[i], "n2": pairs2[i]})
            scores.append(weighted_freq)

    scores = np.array(scores)

    n = len(nodes)
    edges = np.zeros((n, n))
    costs = np.zeros((n, n))

    for i in range(n):
        for j in range(i):
            if nodes[i]["n1"] == nodes[j]["n1"] or \
               nodes[i]["n2"] == nodes[j]["n2"]:
                continue

            is_connected = False

            if p1.edges[nodes[i]["n1"], nodes[j]["n1"]] > 0:
                is_connected = True
            d1 = dist1[nodes[i]["n1"], nodes[j]["n1"]]

            if p2.edges[nodes[i]["n2"], nodes[j]["n2"]] > 0:
                is_connected = True
            d2 = dist2[nodes[i]["n2"], nodes[j]["n2"]]

            if math.fabs(d1 - d2) <= dist_tol:
                if is_connected:
                    costs[i, j] = costs[j, i] = math.fabs(d1 - d2)
                edges[i, j] = edges[j, i] = 1.0

    alignment = list(range(old_len))
    score = np.sum(scores[alignment])
    cost = np.sum(costs[alignment, :][:, alignment]) / 2
    scorecost = score-cost

    for clique in __BronKerbosch(edges):
        clique = list(clique)
        s = np.sum(scores[clique])
        c = np.sum(costs[clique, :][:, clique]) / 2

        if (s - c > scorecost) or (s - c == scorecost and s > score):
            idx1 = []
            idx2 = []
            n1 = []
            n2 = []
            score = s
            cost = c
            scorecost = s - c

            for pair in clique:
                idx1.append(nodes[pair]["n1"])
                idx2.append(nodes[pair]["n2"])
    return score, cost, [idx1, idx2]


def __add_neighbours(p1, p2, n1, n2, idx1, idx2, mapping=None, dist1=None,
                     dist2=None, dist_tol=0, pairs=None):
    """Try to extend alignment by adding nodes connected to already aligned
    parts of pharmacophores. In most cases there is nothing to add, but
    sometimes differences beteween scaffolds of the molecules (rings vs
    linear fragments) result in incomplete alignment after rings decompression.
    Also, pairs of neighbours incompatible with global constraints used
    in coarse-grained alignment (to speed-up the procedure) will be added here.


    Args:
       p1, p2 (Pharmacophore): models to align
       n1, n2 (list of ints): nodes to align
       idx1, idx2 (list of ints): lists of aligned nodes
       mapping (numpy array, optional): array describing nodes compatibility

    Returns:
       float: unnormalized similarity score
       float: edge length differences cost
       2D list: list of two lists representing matched nodes
    """

    if mapping is None:
        mapping = np.zeros((p1.numnodes, p2.numnodes))

        for i in range(p1.numnodes):
            for j in range(p2.numnodes):
                weighted_freq, _ = compare_nodes(p1.nodes[i], p2.nodes[j])
                if weighted_freq > 0.0:
                    mapping[i][j] = weighted_freq
    if dist1 is None:
        dist1 = distances(p1)
        dist1[p1.edges > 0] = p1.edges[p1.edges > 0]

    if dist2 is None:
        dist2 = distances(p2)
        dist2[p2.edges > 0] = p2.edges[p2.edges > 0]

    if pairs is None:
        pairs = []

    def is_compatible(pair1, pair2):
        if (pair1[0] == pair2[0]) or (pair1[1] == pair1[1]):
            return False
        if (p1.edges[pair1[0], pair2[0]] != 0) and \
           (p2.edges[pair1[1], pair2[1]] != 0):
            dist_diff = math.fabs((p1.edges[pair1[0], pair2[0]] -
                                   p2.edges[pair1[1], pair2[1]]))
            if dist_diff <= dist_tol:
                return True
            else:
                return False
        else:
            return True

    for i1, i2 in zip(idx1, idx2):
        neighbours1 = []
        for node in reversed(np.where(p1.edges[i1, n1] > 0)[0]):
            neighbours1.append(n1[node])
            n1.remove(n1[node])

        neighbours2 = []
        for node in reversed(np.where(p2.edges[i2, n2] > 0)[0]):
            neighbours2.append(n2[node])
            n2.remove(n2[node])

        # find compatible neighbours
        for neigh1 in neighbours1:
            for neigh2 in neighbours2:
                if (mapping[neigh1, neigh2] > 0):
                    dist_diff = np.abs(dist1[idx1, neigh1] -
                                       dist2[idx2, neigh2])

                    # there is edge between nodes in at least one of phars
                    connected = np.where(p1.edges[idx1, neigh1] +
                                         p2.edges[idx2, neigh2])
                    max_cost = np.max(dist_diff[connected])
                    if (max_cost <= dist_tol):
                        pairs.append((neigh1, neigh2))

    score = np.sum(mapping[idx1, idx2])
    dist_diff = np.abs(dist1[idx1][:, idx1] - dist2[idx2][:, idx2])
    connected = np.where(p1.edges[idx1][:, idx1] + p2.edges[idx2][:, idx2])

    if len(connected[0]) > 0:
        assert np.max(dist_diff[connected]) <= dist_tol, "cost too high"
    cost = np.sum(dist_diff[connected]) / 2.0

    aln = [idx1[:], idx2[:]]

    for i, pair in enumerate(pairs):
        compatible = [p for p in pairs[i+1:] if is_compatible(pair, p)]

        s, c, (aln1, aln2) = __add_neighbours(p1, p2, n1[:], n2[:],
                                              idx1+[pair[0]], idx2+[pair[1]],
                                              mapping, dist1, dist2, dist_tol,
                                              compatible)

        if (s - c > score - cost) or (s - c == score - cost and s > score):
            score = s
            cost = c
            aln = [aln1[:], aln2[:]]

    return score, cost, aln


def map_pharmacophores(p1, p2, dist_tol=0.0, coarse_grained=True,
                       add_neighbours=False):
    """Find best common substructure match for two Pharmacophores.

    Args:
       p1, p2 (Pharmacophore): models to align
       dist_tol (float, optional): accept distance differences below this
         threshold
       coarse_grained (bool, optional): if True, find alignment for compressed
         ring systems. Otherwise align ring members afer finding coarse-grained
         alignment and try to add neighbours to already aligned nodes.
       add_neighbours (bool, optional): if True, try to extend fine-grained
         alignment by adding neighbours of already aligned nodes. This option
         is ignored if coarse_grained is set to True.

    Returns:
       float: unnormalized similarity score
       float: edge length differences cost
       2D list: list of two lists representing matched nodes. If coarse-grained
         alignment is computed, ring systems are represented as numpy arrays
         containing indices of nodes forming a system
    """
    if not isinstance(p1, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p1).__name__)

    if not isinstance(p2, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p2).__name__)

    if not isinstance(dist_tol, (int, float)):
        raise TypeError("dist_tol must be float or int!")

    if dist_tol < 0:
        raise ValueError("dist_tol must be greater than or equal 0")

    if not isinstance(coarse_grained, bool):
        raise TypeError("coarse_grained must be bool!")

    if not isinstance(add_neighbours, bool):
        raise TypeError("add_neighbours must be bool!")

    if coarse_grained and add_neighbours:
        warnings.warn("Neighbours cannot be added to coarse-grained alignment."
                      "If you want to add neighbours, use coarse_grained=False")

    mapping = np.zeros((p1.numnodes, p2.numnodes))

    for i in range(p1.numnodes):
        for j in range(p2.numnodes):
            weighted_freq, _ = compare_nodes(p1.nodes[i], p2.nodes[j])
            if weighted_freq > 0.0:
                mapping[i][j] = weighted_freq

    dist1 = distances(p1)
    dist1[p1.edges > 0] = p1.edges[p1.edges > 0]
    dist2 = distances(p2)
    dist2[p2.edges > 0] = p2.edges[p2.edges > 0]

    idx1 = [[]]
    idx2 = [[]]
    n1 = [[]]
    n2 = [[]]
    score = 0.0
    cost = 0.0

    scorecost = float("-inf")

    nodes, scores, edges, costs = __modular_product(p1, p2, dist1, dist2,
                                                    dist_tol)

    ring_pairs = [i for i in range(len(nodes)) if "members" in nodes[i]]

    for clique in __BronKerbosch(edges):
        clique = list(clique)
        s = np.sum(scores[clique])
        c = np.sum(costs[clique, :][:, clique]) / 2

        if (s - c >= scorecost):
            if (s - c > scorecost) or (s - c == scorecost and s > score):
                # replace current solutions with a better one
                score = s
                cost = c
                scorecost = s - c

                idx1 = [[]]
                idx2 = [[]]
                n1 = [[]]
                n2 = [[]]

            else:
                # remember another solution with same score
                idx1.append([])
                idx2.append([])
                n1.append([])
                n2.append([])

            rings = [i for i in clique if i in ring_pairs]

            for pair in clique:
                if pair in rings:
                    n1[-1].append(np.array(nodes[pair]["members"][0]))
                    n2[-1].append(np.array(nodes[pair]["members"][1]))
                else:
                    idx1[-1].append(nodes[pair]["n1"])
                    idx2[-1].append(nodes[pair]["n2"])

    if not coarse_grained:
        score = 0.0
        cost = 0.0
        scorecost = float("-inf")
        aln1 = []
        aln2 = []

        for i in range(len(idx1)):
            s, c, [tmp1, tmp2] = __align_rings(p1, p2, n1[i], n2[i],
                                               idx1[i], idx2[i],
                                               mapping, dist1, dist2,
                                               dist_tol)
            if (s - c > scorecost) or (s - c == scorecost and s > score):
                score = s
                cost = c
                scorecost = s - c
                aln1 = tmp1[:]
                aln2 = tmp2[:]

            if add_neighbours:
                remaining1 = [node for node in range(p1.numnodes)
                              if node not in aln1]

                remaining2 = [node for node in range(p2.numnodes)
                              if node not in aln2]

                score, cost, (aln1, aln2) = __add_neighbours(p1, p2,
                                                             remaining1,
                                                             remaining2,
                                                             aln1[:], aln2[:],
                                                             mapping,
                                                             dist1, dist2,
                                                             dist_tol)

    else:
        aln1 = idx1[0]
        aln2 = idx2[0]
        aln1 += n1[0]
        aln2 += n2[0]

    return score, cost, [aln1, aln2]


def similarity(p1, p2, dist_tol=0.0, coarse_grained=True, add_neighbours=False):
    """Find common part of two Pharmacophores, calculate normalized similarity
    score and edge length differences cost.

    Args:
       p1, p2 (pharmacophore): models to align
       dist_tol (float, optional): accept distance differences below this
         threshold
       coarse_grained (bool, optional): if True, find alignment for compressed
         ring systems. Otherwise align ring members afer finding coarse-grained
         alignment and try to add neighbours to already aligned nodes.
       add_neighbours (bool, optional): if True, try to extend fine-grained
         alignment by adding neighbours of already aligned nodes. This option
         is ignored if coarse_grained is set to True.


    Returns:
        score (float): normalized similarity score (value between 0 and 1)
        cost (float): edge length differences cost
    """
    if not isinstance(p1, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p1).__name__)

    if not isinstance(p2, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p2).__name__)

    if not isinstance(dist_tol, (int, float)):
        raise TypeError("dist_tol must be float or int!")

    if dist_tol < 0:
        raise ValueError("dist_tol must be greater than or equal 0")

    if not isinstance(coarse_grained, bool):
        raise TypeError("coarse_grained must be bool!")

    if not isinstance(add_neighbours, bool):
        raise TypeError("add_neighbours must be bool!")

    if p1.numnodes == 0 and p2.numnodes == 0:
        warnings.warn("Pharmacophores are empty!")

    score, cost, _ = map_pharmacophores(p1, p2, dist_tol, coarse_grained,
                                        add_neighbours)
    a1 = 0.0
    a2 = 0.0
    for n in p1.nodes:
        a1 += n["freq"]
    for n in p2.nodes:
        a2 += n["freq"]
    return (score / (a1 + a2)), cost


def combine_pharmacophores(p1, p2, dist_tol=0.0, freq_cutoff=0.0,
                           add_neighbours=False):
    """Create new model from Pharmacophores p1 and p2

    Find common part of two Pharmacophores, add unique elements and calculate
      new frequencies and distances.

    Args:
      p1, p2 (Pharmacophore): models to combine
      dist_tol (float, optional): accept distance differences below this
        threshold
      freq_cutoff (float, optional): skip unique nodes with frequencies below
        this threshold
      add_neighbours (bool, optional): if True, try to extend alignment by
        adding neighbours of already aligned nodes.

    Returns:
       Pharmacophore: combination of p1 and p2
    """
    if not isinstance(p1, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p1).__name__)

    if not isinstance(p2, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p2).__name__)

    if not isinstance(dist_tol, (int, float)):
        raise TypeError("dist_tol must be float or int!")

    if dist_tol < 0:
        raise ValueError("dist_tol must be greater than or equal 0")

    if not isinstance(freq_cutoff, (int, float)):
        raise TypeError("freq_cutoff must be float or int!")

    if freq_cutoff < 0 or freq_cutoff > 1:
        raise ValueError("Invalid freq_cutoff! Use value in the range [0,1]")

    if not isinstance(add_neighbours, bool):
        raise TypeError("add_neighbours must be bool!")

    # find common pharmacophore
    _, _, mapped_nodes = map_pharmacophores(p1, p2, dist_tol,
                                            coarse_grained=False,
                                            add_neighbours=add_neighbours)
    dist1 = distances(p1)
    dist1[p1.edges > 0] = p1.edges[p1.edges > 0]
    dist2 = distances(p2)
    dist2[p2.edges > 0] = p2.edges[p2.edges > 0]

    # we will need it later
    added = {0: {}, 1: {}}

    # create new graph from common part
    molecules = p1.molecules + p2.molecules

    title = "("+p1.title+")+("+p2.title+")"
    nodes = []

    idx = 0
    for i in range(len(mapped_nodes[0])):
        u = p1.nodes[mapped_nodes[0][i]]
        v = p2.nodes[mapped_nodes[1][i]]
        _, types = compare_nodes(u, v)
        nodes.append({"label": idx, "type": types,
                      "freq": u["freq"] + v["freq"]})
        for j in [0, 1]:
            added[j][idx] = mapped_nodes[j][i]
        idx += 1

    # add edges
    edges = np.zeros((idx, idx))
    for i in range(idx):
        no1 = (added[0][i], added[1][i])
        for j in range(i):
            dist = 0.0
            no2 = (added[0][j], added[1][j])
            freq1 = p1.nodes[no1[0]]["freq"] + p1.nodes[no2[0]]["freq"]
            freq2 = p2.nodes[no1[1]]["freq"] + p2.nodes[no2[1]]["freq"]
            if p1.edges[no1[0], no2[0]] or p2.edges[no1[1], no2[1]]:
                d1 = dist1[no1[0], no2[0]]
                d2 = dist2[no1[1], no2[1]]
                dist = (d1 * freq1 + d2 * freq2) / (freq1 + freq2)
                edges[i, j] = edges[j, i] = dist

    # do not warn about empty pharmacophore (nodes might be added latter)
    # warn about empty common part instead
    warnings.simplefilter("ignore", UserWarning)
    new_p = Pharmacophore(nodes=nodes, edges=edges, molecules=molecules,
                          title=title)
    warnings.simplefilter("always", UserWarning)
    if new_p.numnodes == 0:
        warnings.warn("Empty common part!")

    # add unique elements
    freq_cutoff = molecules * freq_cutoff

    to_add = [[i for i in range(p1.numnodes) if i not in mapped_nodes[0] and
               p1.nodes[i]["freq"] >= freq_cutoff],
              [i for i in range(p2.numnodes) if i not in mapped_nodes[1] and
               p2.nodes[i]["freq"] >= freq_cutoff]]

    for (nr, phar) in {0: p1, 1: p2}.items():
        for n in to_add[nr]:
            added[nr][idx] = n
            new_p.add_node(phar.nodes[n].copy())
            new_p.nodes[idx]["label"] = idx
            for (k, v) in added[nr].items():
                if phar.edges[n, v] > 0:
                    new_p.add_edge(k, idx, phar.edges[n, v])
            idx += 1

    # check if new pharmacophore is connected
    components = split_components(new_p)
    comp_nr = len(components)

    if comp_nr > 1:
        # shortest distances between components
        comp_dist = np.zeros((comp_nr, comp_nr)) + float("inf")

        # nearest_node[i, j] == id of node from component j, that is nearest to
        # component i
        nearest_node = np.zeros((comp_nr, comp_nr), dtype=int)
        for i in range(comp_nr):
            for j in range(i):
                shortest_dist = float("inf")
                nearest_nodes = [None, None]
                for n1 in components[i]:
                    for n2 in components[j]:
                        if n1 in added[0] and n2 in added[0]:
                            d1 = dist1[added[0][n1], added[0][n2]]
                            freq1 = p1.nodes[added[0][n1]]["freq"] + \
                                    p1.nodes[added[0][n2]]["freq"]
                        else:
                            d1 = 0
                            freq1 = 0
                        if n1 in added[1] and n2 in added[1]:
                            d2 = dist2[added[1][n1], added[1][n2]]
                            freq2 = p2.nodes[added[1][n1]]["freq"] + \
                                    p2.nodes[added[1][n2]]["freq"]
                        else:
                            d2 = 0
                            freq2 = 0
                        if (freq1 + freq2) == 0:
                            dist = float("inf")
                        else:
                            dist = (d1 * freq1 + d2 * freq2) / (freq1 + freq2)

                        if dist < shortest_dist:
                            shortest_dist = dist
                            nearest_nodes = [n1, n2]
                comp_dist[i, j] = comp_dist[j, i] = shortest_dist
                if shortest_dist < float("inf"):
                    nearest_node[i, j] = nearest_nodes[1]
                    nearest_node[j, i] = nearest_nodes[0]

        sorted_connections = np.unravel_index(comp_dist.argsort(axis=None),
                                              comp_dist.shape)

        # connect components
        for i, j in zip(*sorted_connections):
            n1 = nearest_node[i, j]
            n2 = nearest_node[j, i]
            new_p.add_edge(n1, n2, comp_dist[i, j])

            # check if graph is already connected
            if len(split_components(new_p)) == 1:
                break

    return new_p


def inclusive_similarity(p1, p2, dist_tol=0.0, coarse_grained=True,
                         add_neighbours=False):
    """Find common part of two Pharmacophores and calculate what fractions of
    both models it contains. E.g. if p1 is a substructure of p2, function will
    return (1.0, s, c), where s<=1 and c>=0.

    Args:
       p1, p2 (pharmacophore): models to align
       dist_tol (float, optional): accept distance differences below this
         threshold
       coarse_grained (bool, optional): if True, find alignment for compressed
         ring systems. Otherwise align ring members afer finding coarse-grained
         alignment and try to add neighbours to already aligned nodes.
       add_neighbours (bool, optional): if True, try to extend fine-grained
         alignment by adding neighbours of already aligned nodes. This option
         is ignored if coarse_grained is set to True.


    Returns:
        float: normalized similarity score for p1
        float: normalized similarity score for p2
        float: edge length differences cost
    """
    if not isinstance(p1, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p1).__name__)

    if not isinstance(p2, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p2).__name__)

    if not isinstance(dist_tol, (int, float)):
        raise TypeError("dist_tol must be float or int!")

    if dist_tol < 0:
        raise ValueError("dist_tol must be greater than or equal 0")

    if not isinstance(coarse_grained, bool):
        raise TypeError("coarse_grained must be bool!")

    if not isinstance(add_neighbours, bool):
        raise TypeError("add_neighbours must be bool!")

    score, cost, _ = map_pharmacophores(p1, p2, dist_tol, coarse_grained,
                                        add_neighbours)
    a1 = 0.0
    a2 = 0.0
    m1 = p1.molecules
    m2 = p2.molecules
    m = m1 + m2
    for n in p1.nodes:
        a1 += n["freq"]
    for n in p2.nodes:
        a2 += n["freq"]
    return score * (m1 / m) / a1, score * (m2 / m) / a2, cost


def filter_nodes(p, freq_range=(0.0, 1.0), rm_outside=True):
    """Create new model without nodes that does not fulfill given frequency
    criteria.

    Args:
       p (Pharmacophore): model to filter
       freq_range (tuple, optional): two floats, frequence range for filtering
       rm_outside (bool, optional): if True remove nodes with frequencies
         outside given range; remove nodes with frequencies inside range
         otherwise

    Returns:
       Pharmacohpre: new model
    """

    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    if not isinstance(freq_range, tuple):
        raise TypeError("Invalid freq_range!")

    elif not ((len(freq_range) == 2)
              and isinstance(freq_range[0], (float, int))
              and isinstance(freq_range[1], (float, int))
              and (freq_range[0] <= freq_range[1])):
        raise ValueError("Invalid freq_range!")
    elif not (freq_range[0] >= 0 and freq_range[1] <= 1):
        raise ValueError("Invalid freq_range! Use values in the range [0,1]")

    if not isinstance(rm_outside, bool):
        raise TypeError("rm_outside should be bool!")

    dist = distances(p)
    dist[p.edges > 0] = p.edges[p.edges > 0]
    new_p = p.copy()

    freq_range = [i * p.molecules for i in freq_range]
    if rm_outside:
        for n in range(p.numnodes - 1, -1, -1):
            if new_p.nodes[n]["freq"] < freq_range[0] or \
               new_p.nodes[n]["freq"] > freq_range[1]:
                new_p.remove_node(n)
                dist = np.delete(np.delete(dist, n, 0), n, 1)

    else:
        for n in range(p.numnodes - 1, -1, -1):
            if new_p.nodes[n]["freq"] >= freq_range[0] or \
               new_p.nodes[n]["freq"] <= freq_range[1]:
                new_p.remove_node(n)
                dist = np.delete(np.delete(dist, n, 0), n, 1)

    # check if new pharmacophore is connected
    components = split_components(new_p)
    comp_nr = len(components)

    if comp_nr > 1:
        # shortest distances between components
        comp_dist = np.zeros((comp_nr, comp_nr)) + float("inf")

        # nearest_node[i, j] == id of node from component j, that is nearest to
        # component i
        nearest_node = np.zeros((comp_nr, comp_nr), dtype=int)
        for i in range(comp_nr):
            for j in range(i):
                shortest_dist = float("inf")
                nearest_nodes = [None, None]
                for n1 in components[i]:
                    for n2 in components[j]:
                        if dist[n1, n2] < shortest_dist:
                            shortest_dist = dist[n1, n2]
                            nearest_nodes = [n1, n2]
                comp_dist[i, j] = comp_dist[j, i] = shortest_dist
                nearest_node[i, j] = nearest_nodes[1]
                nearest_node[j, i] = nearest_nodes[0]

        shortest_connection = np.argmin(comp_dist, axis=1)

        # connect components
        for i in range(comp_nr):
            j = shortest_connection[i]
            n1 = nearest_node[i, j]
            n2 = nearest_node[j, i]
            new_p.add_edge(n1, n2, comp_dist[i, j])

    return new_p


def spring_layout(p, c0=0.2, c1=1.0):
    """Calculate points positions for Pharmacophore depiction using spring
    layout.

    Args:
       p (Pharmacophore): model to depict
       c0, c1 (float, optional): coefficients for spring and repulsive forces

    Returns:
       numpy array: 2D array with nodes positions
    """
    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    if not (isinstance(c0, (float, int)) and isinstance(c1, (float, int))):
        raise TypeError("Invalid constants!")
    if not (c0 > 0 and c1 > 0):
        raise ValueError("Invalid constants! Use values greater than 0.")

    if p.numnodes == 0:
        raise ValueError("Pharmacophore is empty!")

    from scipy.optimize import minimize
    from scipy.spatial.distance import pdist, squareform

    def f(x):
        x = x.reshape(int(len(x)/2), 2)
        eng = 0
        norms = pdist(x)
        norms[norms == 0] = 0.000001
        e = squareform(p.edges)

        spring = (norms - e)[np.where(e > 0)]**2
        eng += np.sum(spring) * c1

        repulsive = norms[np.where(e == 0)]
        eng += np.sum(c0 / repulsive[np.nonzero(repulsive)])
        return eng

    x0 = np.random.random(p.numnodes*2)

    res = minimize(f, x0)
    newpositions = res.x.reshape((p.numnodes, 2))
    return newpositions


def draw(p, layout="rd"):
    """Draw Pharmacophore using RDKit ("rd"), OpenBabel ("ob") or spring
    layout ("spring") to calculate nodes positions.

    We recommend to use RDKit ("rd"), as it generates the clearest layouts.

    Args:
       p (Pharmacophore): model to depict
       layout (str, optional): layout name

    Returns:
       tuple:
         * matplotlib Figure
         * matplotlib axis
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge
    from matplotlib.font_manager import FontManager

    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore, got %s instead" %
                        type(p).__name__)

    if not isinstance(layout, str):
        raise TypeError("Invalid layout! Expected str, got %s instead." %
                        type(layout).__name__)

    if p.numnodes == 0:
        raise ValueError("Pharmacophore is empty!")

    if layout == "rd":
        try:
            from decaf.toolkits.rd import layout
            pos = layout(p)
        except Exception as e:
            raise ImportError("Cannot use 'rd' layout! Use 'ob' or 'spring'"
                              "instead", e)

    elif layout == "ob":
        try:
            from decaf.toolkits.ob import layout
            pos = layout(p)
        except Exception as e:
            raise ImportError("Cannot use 'ob' layout! Use 'rd' or 'spring'"
                              "instead", e)

    elif layout == "spring":
        try:
            pos = spring_layout(p)
        except Exception as e:
            raise ImportError("Cannot use spring layout!", e)
    else:
        raise ValueError("Wrong layout specified! Use 'rd', 'ob' or 'spring'"
                         "instead.")

    ax_coeff = 1.

    def fontsize(idx, default=FontManager.get_default_size()):
        coeff = p.nodes[idx]["freq"] / p.molecules
        size = default * coeff * ax_coeff
        return size

    fig, ax = plt.subplots()
    plt.axis("equal")
    plt.axis("off")

    axis = (np.min(pos[:, 0])-1,
            np.max(pos[:, 0])+1,
            np.min(pos[:, 1])-1,
            np.max(pos[:, 1])+1)
    plt.axis(axis)

    # calculate scaling ratio for font
    ax_coeff = 12. / max((axis[1]-axis[0]), (axis[3]-axis[2]))

    for i in range(p.numnodes):
        for j in range(i):
            if p.edges[i, j] > 0:
                tmp = np.array([pos[i], pos[j]])
                ax.plot(tmp[:, 0], tmp[:, 1], color="#000000", zorder=1)

        r = p.nodes[i]["freq"] / p.molecules * 0.3
        fsize = fontsize(i)
        nfreq = sum(p.nodes[i]["type"].values())
        theta1 = 0.0
        for t in p.nodes[i]["type"]:
            delta = 360 * p.nodes[i]["type"][t] / nfreq
            theta2 = theta1+delta
            w = Wedge(pos[i], r, theta1, theta2, ec="none", fc=COLORS[t])
            ax.add_artist(w)
            ax.text(pos[i][0], pos[i][1], str(p.nodes[i]["label"]),
                    color="#000000", ha="center", va="center", size=fsize)
            theta1 = theta2

    plt.show()
    return fig, ax
