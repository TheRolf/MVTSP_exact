import itertools
import operator

from mvtsp_exact.tools.degree_sequences import combination_to_sequence


def nonNegativeQ(vector):
    return all(i >= 0 for i in vector)


def tree_U_multigraph(x_st, x_tp):
    for i, j in x_st:
        if (i, j) not in x_tp:
            x_tp[i, j] = 0
        x_tp[i, j] += 1
    return x_tp


def distribute_degree(Do, Di, k, o_max, i_max):
    if k==1:
        if Do <= o_max[0] and Di <= i_max[0]:
            yield [Do], [Di]
    else:
        if Do==0 and Di>0:
            ds_o = [0]*k
            for comb2 in itertools.combinations(list(range(k-1+Di)), k-1):
                ds_i = combination_to_sequence(comb2, Di)
                if nonNegativeQ(list(map(operator.sub, i_max, ds_i))):
                    yield ds_o, ds_i

        if Do>0 and Di==0:
            ds_i = [0]*k
            for comb1 in itertools.combinations(list(range(k-1+Do)), k-1):
                ds_o = combination_to_sequence(comb1, Do)
                if nonNegativeQ(list(map(operator.sub, o_max, ds_o))):
                    yield ds_o, ds_i

        if Do>0 and Di>0:
            for comb1 in itertools.combinations(list(range(k-1+Do)), k-1):
                ds_o = combination_to_sequence(comb1, Do)
                for comb2 in itertools.combinations(list(range(k-1+Di)), k-1):
                    ds_i = combination_to_sequence(comb2, Di)
                    if nonNegativeQ(list(map(operator.sub, o_max, ds_o))) and \
                       nonNegativeQ(list(map(operator.sub, i_max, ds_i))):
                        yield ds_o, ds_i


def subsets(baseset, emptyQ=False, minSize=None, maxSize=None, fullQ=False):
    n = len(baseset)
    minSize = (0 if emptyQ else 1) if minSize is None else minSize
    maxSize = (n if fullQ else n-1) if maxSize is None else maxSize
    for size in range(minSize, maxSize+1):
        for s in itertools.combinations(baseset, size):
            yield set(s)


def subset_pairs(baseset, emptyFullQ=False, minSize=None, maxSize=None):
    if minSize is not None and maxSize is not None and (2*minSize > len(baseset) or 2*maxSize < len(baseset)):
        yield [], []
    else:
        n = len(baseset)
        baseset = set(baseset)
        minSize = (0 if emptyFullQ else 1) if minSize is None else minSize
        maxSize = (n if emptyFullQ else n-1) if maxSize is None else maxSize
        for size in range(minSize, maxSize + 1):
            for s1 in itertools.combinations(baseset, size):
                s2 = baseset.difference(s1)
                yield (set(s1), set(s2))
