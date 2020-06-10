import operator

from mvtsp_exact.tools.degree_sequences import directed_tree_DS_wrapper, directed_DS_to_tree
from mvtsp_exact.tools.tools import nonNegativeQ, tree_U_multigraph
from mvtsp_exact.tools.transport import transport


def enum_DS(c, r, rootedQ=False, boundedQ=False):
    """Solves the MV-TSP problem using enumerating directed trees (improved)

    Input: cost matrix 'c', requiremets vector 'r'
    Output: objective value, dictionary of edges with multiplicities

    The algorithm generates all c^N degree sequences of directed trees on N vertices, and solves the underlying
    transportation problem. For each degree sequence, it enumerates all directed spanning trees realising the degree
    sequence, and chooses the cheapest one. It returns the multigraph with minimal overall cost.
    Three modes:
    1) rootedQ=False, boundedQ=False: it enumerates degree sequences of directed trees (~11.09^n)
    2) rootedQ=True, boundedQ=False: it enumerates degree sequences of rooted trees (~4^n)
    3) rootedQ=True, boundedQ=True: it enumerates degree sequences of rooted trees, but omits infeasible ones (i.e.
       where d_T(v) would be higher than r(v)
    Note: in the paper, this algorithm is described in §2.2.1 as an improved version of enum-MV

    Time complexity:  n^n (improved)
    Space complexity: poly(n)
    """
    N = len(r)
    best_obj = float('inf')
    best_tp = []
    best_st = []
    bounds = [r[i]-1 if i==0 else r[i] for i in range(N)] if boundedQ else None
    for sb1, sb2 in directed_tree_DS_wrapper(N, rootedQ=rootedQ, bounds=bounds):
        d_out = list(map(operator.sub, r, sb1))
        d_in = list(map(operator.sub, r, sb2))
        if nonNegativeQ(d_out) and nonNegativeQ(d_in):
            obj_tp, x_tp = transport(c, d_out, d_in)
            obj_st, x_st = float('inf'), None
            for T in directed_DS_to_tree(sb1, sb2, []):
                obj_st_loc = 0
                for u, v in T:
                    obj_st_loc += c[u][v]
                if obj_st_loc < obj_st:
                    obj_st = obj_st_loc
                    x_st = T
            if obj_st + obj_tp < best_obj:
                best_obj = obj_st + obj_tp
                best_tp = x_tp
                best_st = x_st
    best_sol = tree_U_multigraph(best_st, best_tp)
    return best_obj, best_sol


if __name__ == '__main__':
    pass
