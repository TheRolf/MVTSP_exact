from math import floor, ceil, log
from operator import sub, add

# from src.math_tools import nonNegativeQ
# from src.mvtsp.mvtsp_tools import subSets, subSetPairs, multigraph
# from src.mvtsp.spanning_tree_tools import directedTreeDSwrapper, degreeDist
# from src.mvtsp.transport import transport
# from src.tools import dict_union
from mvtsp_exact.tools.degree_sequences import directed_tree_DS_wrapper
from mvtsp_exact.tools.tools import nonNegativeQ, tree_U_multigraph, distribute_degree, subset_pairs, subsets
from mvtsp_exact.tools.transport import transport


def dc_MV(c, r, perfectlyBalancedQ=False, rootedQ=False, boundedQ=False):
    """Solves the MV-TSP problem using divide & conquer

    Input: cost matrix 'c', requiremets vector 'r'
    Output: objective value, dictionary of edges with multiplicities

    The algorithm generates all c^N degree sequences of directed trees on N vertices, and solves the underlying
    transportation problem. For each degree sequence, it calculates the minimum cost spanning tree realising that
    degree sequence, using a divide & conquer approach inspired by Gurevich & Shelah.
    Options:
    - rootedQ: considers rooted spanning trees, if set True
    - boundedQ: does not generate infeasible degree sequences, if set True (only works if rootedQ=True)
    - perfectlyBalancedQ: uses perfectly balanced (~1/2-1/2) partitioning, if set True; ~1/3-2/3 otherwise
      (only works if rootedQ=True)

    Complexities if both rootedQ and perfectlyBalancedQ are set True
    Time complexity:  16^n
    Space complexity: poly(n)
    """
    N = len(r)
    best_obj = float('inf')
    best_tp = []
    best_st = []
    cdict = {(i, j): c[i][j] for i in range(N) for j in range(N) if i != j}
    bounds = [r[i]-1 if i==0 else r[i] for i in range(N)] if boundedQ else None
    for sb1, sb2 in directed_tree_DS_wrapper(N, rootedQ=rootedQ, bounds=bounds):
        sb1dict, sb2dict = {i: sb1[i] for i in range(N)}, {i: sb2[i] for i in range(N)}
        d_out, d_in = list(map(sub, r, sb1)), list(map(sub, r, sb2))
        if nonNegativeQ(d_out) and nonNegativeQ(d_in):
            obj_tp, x_tp = transport(c, d_out, d_in)
            obj_st, x_st = T(cdict, set(range(N)), sb1dict, sb2dict, perfectlyBalancedQ=perfectlyBalancedQ)
            if obj_st + obj_tp < best_obj:
                best_obj = obj_st + obj_tp
                best_tp = x_tp
                best_st = x_st
    best_sol = tree_U_multigraph(best_st, best_tp)
    return best_obj, best_sol


def T(c, S, d_o, d_i, perfectlyBalancedQ=False):
    if len(S) == 2:
        return cost_of_two(d_o, d_i, c)

    elif len(S) == 3:
        return cost_of_three(d_o, d_i, c)

    else:
        cost = float("inf")
        tree = []
        s0 = min(S) - 1

        if len(S) <= 9:
            perfectlyBalancedQ = False

        minSize = floor(len(S)/2) if perfectlyBalancedQ else floor(len(S)/3)
        maxSize = ceil(len(S)/2)if perfectlyBalancedQ else ceil(2*len(S)/3)

        if perfectlyBalancedQ:
            for S1, S2 in subset_pairs(S, minSize=minSize, maxSize=maxSize):
                for seps in subsets(S1, minSize=1, maxSize=floor(log(len(S)))):
                    k = len(seps)
                    d_o_sep_in_S1 = len(S1)-1 - sum(d_o[p] for p in S1 if p not in seps)
                    d_i_sep_in_S1 = len(S1)-1 - sum(d_i[p] for p in S1 if p not in seps)
                    d_o_sep_in_S2 = len(S2)+k-1 - sum(d_o[p] for p in S2)
                    d_i_sep_in_S2 = len(S2)+k-1 - sum(d_i[p] for p in S2)
                    ds_o_sep, ds_i_sep = [d_o[p] for p in seps], [d_i[p] for p in seps]

                    if nonNegativeQ([d_o_sep_in_S1, d_i_sep_in_S1, d_o_sep_in_S2, d_i_sep_in_S2]):
                        for ds_o_sep_in_S1, ds_i_sep_in_S1 in distribute_degree(
                                d_o_sep_in_S1, d_i_sep_in_S1, k, ds_o_sep, ds_i_sep):
                            ds_o_sep_in_S2 = list(map(sub, ds_o_sep, ds_o_sep_in_S1))
                            ds_i_sep_in_S2 = list(map(sub, ds_i_sep, ds_i_sep_in_S1))

                            if ds_i_sep_in_S2.count(1) <= 1:
                                if ds_i_sep_in_S2.count(1) == 0:
                                    ds_o_s0_in_S2 = k
                                    ds_i_s0_in_S2 = 0
                                    for i in range(len(seps)):
                                        ds_i_sep_in_S2[i] += 1

                                else:  # ds_i_sep_in_S2.count(1) == 1
                                    r = ds_i_sep_in_S2.index(1)
                                    ds_o_sep_in_S2[r] += 1
                                    for i in range(len(seps)):
                                        if i!=r:
                                            ds_i_sep_in_S2[i] += 1
                                    ds_o_s0_in_S2 = k-1
                                    ds_i_s0_in_S2 = 1

                                if nonNegativeQ(ds_o_sep_in_S2) and nonNegativeQ(ds_i_sep_in_S2) and \
                                        list(map(add, ds_o_sep_in_S1, ds_i_sep_in_S1)).count(0)==0 and \
                                        list(map(add, ds_o_sep_in_S2, ds_i_sep_in_S2)).count(0)==0:
                                    d_o_1, d_i_1, d_o_2, d_i_2 = {}, {}, {}, {}
                                    d_o_1.update({p: d_o[p] for p in S1 if p not in seps})
                                    d_o_1.update({seps[i]: ds_o_sep_in_S1[i] for i in range(k)})

                                    d_i_1.update({p: d_i[p] for p in S1 if p not in seps})
                                    d_i_1.update({seps[i]: ds_i_sep_in_S1[i] for i in range(k)})

                                    d_o_2.update({p: d_o[p] for p in S2})
                                    d_o_2.update({seps[i]: ds_o_sep_in_S2[i] for i in range(k)})
                                    d_o_2.update({s0: ds_o_s0_in_S2})

                                    d_i_2.update({p: d_i[p] for p in S2})
                                    d_i_2.update({seps[i]: ds_i_sep_in_S2[i] for i in range(k)})
                                    d_i_2.update({s0: ds_i_s0_in_S2})

                                    c1 = {(i, j): c[i, j] for i in S1 for j in S1 if i!=j}
                                    S2_ = S2.union(set(seps))
                                    c2 = {(i, j): c[i, j] for i in S2_ for j in S2_ if i!=j}

                                    c2.update({(i, j): float("inf") for i in seps for j in seps if i!=j})
                                    c2.update({(i, s0): float("inf") for i in S2})
                                    c2.update({(s0, j): float("inf") for j in S2})
                                    c2.update({(i, s0): 0 for i in seps})
                                    c2.update({(s0, j): 0 for j in seps})

                                    cost1, tree1 = T(c1, S1, d_o_1, d_i_1, perfectlyBalancedQ=perfectlyBalancedQ)
                                    cost2, tree2 = T(c2, S2_.union({s0}), d_o_2, d_i_2, perfectlyBalancedQ=perfectlyBalancedQ)  #
                                    if cost1 + cost2 < cost:
                                        cost = cost1 + cost2
                                        tree = tree1 + tree2

        else:  # perfectlyBalancedQ=False
            for S1, S2 in subset_pairs(S, minSize=minSize, maxSize=maxSize):
                for s1 in S1:
                    d_o_s1_in_S1 = len(S1)-1 - sum(d_o[p] for p in S1 if p!=s1)
                    d_i_s1_in_S1 = len(S1)-1 - sum(d_i[p] for p in S1 if p!=s1)
                    d_o_s1_in_S2 = len(S2)+1-1 - sum(d_o[p] for p in S2) # +1: the size of S2 will increase by one (s1)
                    d_i_s1_in_S2 = len(S2)+1-1 - sum(d_i[p] for p in S2)

                    if 0 <= d_o_s1_in_S1 and 0 <= d_i_s1_in_S1 and 0 <= d_o_s1_in_S2 and 0 <= d_i_s1_in_S2 and \
                       1 <= d_o_s1_in_S1+d_i_s1_in_S1 and 1 <= d_o_s1_in_S2+d_i_s1_in_S2:
                        d_o_1 = {p: d_o_s1_in_S1 if p==s1 else d_o[p] for p in S1}
                        d_i_1 = {p: d_i_s1_in_S1 if p==s1 else d_i[p] for p in S1}
                        d_o_2 = {p: d_o[p] for p in S2}
                        d_o_2.update({s1: d_o_s1_in_S2})
                        d_i_2 = {p: d_i[p] for p in S2}
                        d_i_2.update({s1: d_i_s1_in_S2})

                        S2_ = S2.union({s1})

                        cost1, tree1 = T(c, S1, d_o_1, d_i_1, perfectlyBalancedQ=perfectlyBalancedQ)
                        cost2, tree2 = T(c, S2_, d_o_2, d_i_2, perfectlyBalancedQ=perfectlyBalancedQ)
                        if cost1 + cost2 < cost:
                            cost = cost1 + cost2
                            tree = tree1 + tree2
        return cost, tree


def cost_of_two(d_o, d_i, c):
    p_o, p_i = None, None
    for p in d_o:
        if d_o[p] == 1:
            p_o = p
        if d_i[p] == 1:
            p_i = p
    return c[p_o, p_i], [[p_o, p_i]]


def cost_of_three(d_o, d_i, c):
    lf1, lf2, cnt = None, None, None
    types = {}
    for p in d_o:
        if d_o[p]+d_i[p] == 2:
            cnt = p
        elif d_o[p] == 1:
            lf1, lf2 = p, lf1
            types["o"] = p
        elif d_i[p] == 1:
            lf1, lf2 = p, lf1
            types["i"] = p
    if "i" in types and "o" in types:
        return c[types["o"], cnt] + c[cnt, types["i"]], [[types["o"], cnt], [cnt, types["i"]]]
    elif "o" in types:
        return c[lf1, cnt] + c[lf2, cnt], [[lf1, cnt], [lf2, cnt]]
    elif "i" in types:
        return c[cnt, lf1] + c[cnt, lf2], [[cnt, lf1], [cnt, lf2]]


if __name__ == '__main__':
    pass
