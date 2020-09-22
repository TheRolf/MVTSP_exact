import operator

from mvtsp_exact.tools.degree_sequences import directed_tree_DS_wrapper
from mvtsp_exact.tools.tools import nonNegativeQ, tree_U_multigraph
from mvtsp_exact.tools.transport import transport


def dp_MV(c, r, rootedQ=True, boundedQ=True):
    """Solves the MV-TSP problem using dynamic programming

    Input: cost matrix 'c', requiremets vector 'r'
    Output: objective value, dictionary of edges with multiplicities

    The algorithm generates all 2^O(n) degree sequences of directed trees on n vertices, and solves the underlying
    transportation problem. For each degree sequence, it calculates the minimum cost spanning tree realising that
    degree sequence. It returns the multigraph with minimal overall cost.
    1) rootedQ=False, boundedQ=False: it enumerates degree sequences of directed trees (~11.09^n)
    2) rootedQ=True, boundedQ=False: it enumerates degree sequences of rooted trees (~4^n)
    3) rootedQ=True, boundedQ=True: it enumerates degree sequences of rooted trees, but omits infeasible ones (i.e.
       where d_T(v) would be higher than r(v), time and space complexity falls back to O(2^n) for the single-visit TSP

    Complexities for options 2)-3) (log sum r factors omitted)
    Time complexity:  O(5^n)
    Space complexity: O(5^n)
    """
    def T(ds_o, ds_i):
        if (ds_o, ds_i) in Table:
            return Table[(ds_o, ds_i)]

        # continuing: if value is not in the Table, it has to be calculated
        t_new = []
        if sum(ds_o) + sum(ds_i) == 2:
            u = ds_o.index(1)
            v = ds_i.index(1)
            cost = c[u][v] if u != v else float("inf")
            Table[ds_o, ds_i] = (cost, [[u, v]])
            return cost, [[u, v]]

        else:
            val = float("inf")
            for w in [i for i in range(N) if ds_o[i] + ds_i[i] == 1]:
                if ds_o[w] == 1:
                    for v in [i for i in range(N) if ds_i[i] + ds_o[i] > 1 and ds_i[i] > 0 and i != w]:
                        ds_ox = tuple(d-1 if k==w else d for k, d in enumerate(ds_o))
                        ds_ix = tuple(d-1 if k==v else d for k, d in enumerate(ds_i))
                        if c[w][v] + T(ds_ox, ds_ix)[0] < val:
                            val = T(ds_ox, ds_ix)[0] + c[w][v]
                            t_new = T(ds_ox, ds_ix)[1] + [[w, v]]
                if ds_i[w] == 1:
                    for u in [i for i in range(N) if ds_i[i] + ds_o[i] > 1 and ds_o[i] > 0 and i != w]:
                        ds_ox = tuple(d-1 if k==u else d for k, d in enumerate(ds_o))
                        ds_ix = tuple(d-1 if k==w else d for k, d in enumerate(ds_i))
                        if c[u][w] + T(ds_ox, ds_ix)[0] < val:
                            val = T(ds_ox, ds_ix)[0] + c[u][w]
                            t_new = T(ds_ox, ds_ix)[1] + [[u, w]]
            Table[(ds_o, ds_i)] = (val, t_new)
            return val, t_new

    N = len(r)
    best_obj = float('inf')
    best_tp = []
    best_st = []
    Table = {}
    bounds = [r[i]-1 if i==0 else r[i] for i in range(N)] if boundedQ else None
    for sb1, sb2 in directed_tree_DS_wrapper(N, rootedQ=rootedQ, bounds=bounds):
        d_out = list(map(operator.sub, r, sb1))
        d_in = list(map(operator.sub, r, sb2))
        if nonNegativeQ(d_in) and nonNegativeQ(d_out):
            obj_tp, x_tp = transport(c, d_out, d_in)
            obj_st, x_st = T(sb1, sb2)
            if obj_st + obj_tp < best_obj:
                best_obj = obj_st + obj_tp
                best_tp = x_tp
                best_st = x_st
    best_sol = tree_U_multigraph(best_st, best_tp)
    return best_obj, best_sol
