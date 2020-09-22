import itertools

from mvtsp_exact.tools.tools import nonNegativeQ, tree_U_multigraph
from mvtsp_exact.tools.transport import transport


def enum_MV(c, r):
    """Solves the MV-TSP problem using enumerating directed trees

    Input: cost matrix 'c', requiremets vector 'r'
    Output: objective value, dictionary of edges with multiplicities

    The algorithm generates all n^{n-1} Pruefer-codes on N vertices, then all 2^{n-1} orientations of the edges.
    For each such directed spanning tree, it solves a transportation problem to ensure the degree requirements.
    It returns the cheapest of such multigraphs.

    Complexities (log sum r factors omitted)
    Time complexity:  O(n^n)
    Space complexity: poly(n)
    """
    N = len(r)
    best_obj = float('inf')
    best_tp = []
    best_st = []
    for st in directedSpanningTrees(N):
        d_out = r[:]
        d_in = r[:]
        for u, v in st:
            d_out[u] -= 1
            d_in[v] -= 1
        # filter out spanning trees that has vertices with too high out- or indegrees
        if nonNegativeQ(d_out) and nonNegativeQ(d_in):
            obj_st = 0
            for u, v in st:
                obj_st += c[u][v]
            obj_tp, x_tp = transport(c, d_out, d_in)
            if obj_st + obj_tp < best_obj:
                best_obj = obj_st + obj_tp
                best_tp = x_tp
                best_st = st
    best_sol = tree_U_multigraph(best_st, best_tp)
    return best_obj, best_sol


def directedSpanningTrees(N):
    # all possible Pruefer-sequences
    for seq in itertools.product(list(range(N)), repeat=N-2):
        # all possible orientations of edges
        for T in directUndirectedTrees(seq):
            yield T


def directUndirectedTrees(seq):
    # Pruefer-sequence to graph
    # https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
    N = len(seq) + 2
    T = [[] for i in range(N-1)]
    deg = [1 for i in range(N)]
    for i in seq:
        deg[i] += 1
    k = 0
    for i in seq:
        for j in range(N):
            if deg[j] == 1:
                T[k] = [i, j]
                k += 1
                deg[i] -= 1
                deg[j] -= 1
                break
    u, v = 0, 0
    for i in range(N):
        if deg[i] == 1:
            if u == 0:
                u = i
            else:
                v = i
                break
    T[-1] = [v, u]

    # all possible directed trees
    for I in range(2**(N-1)):
        yield [(e[0], e[1]) if (I >> i) & 1 else (e[1], e[0]) for i, e in enumerate(T)]
