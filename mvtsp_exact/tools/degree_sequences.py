import itertools


def directed_DS_to_tree(ds1, ds2, subtree):
    ds = [ds1[i] + ds2[i] for i in range(len(ds1))]
    if sum(ds) == 2:
        i = ds1.index(1)
        j = ds2.index(1)
        yield subtree + [[i, j]]
    else:
        i = ds.index(1)
        if ds1[i] == 1:
            for j in range(len(ds)):
                if i != j and ds2[j] > 0 and ds[j] > 1:
                    for tree in directed_DS_to_tree(
                            [d-1 if k==i else d for k, d in enumerate(ds1)],
                            [d-1 if k==j else d for k, d in enumerate(ds2)],
                            [edge for edge in subtree] + [[i, j]]):
                        yield tree
        else:  # ds2[i] == 1  (i == ds2.index(1))
            for j in range(len(ds)):
                if i != j and ds1[j] > 0 and ds[j] > 1:
                    for tree in directed_DS_to_tree(
                            [d-1 if k==j else d for k, d in enumerate(ds1)],
                            [d-1 if k==i else d for k, d in enumerate(ds2)],
                            [edge for edge in subtree] + [[j, i]]):
                        yield tree


def directed_tree_DS_wrapper(N, rootedQ=False, bounds=None):
    if rootedQ:
        for ds1, ds2 in directed_rooted_tree_DS(N, bounds=bounds):
            yield ds1, ds2
    else:
        for ds1, ds2 in directed_tree_DS(N):
            yield ds1, ds2


def directed_tree_DS(n):
    for comb1 in itertools.combinations(list(range(2*n - 2)), n-1):
        sb1 = combination_to_sequence(comb1, n-1)
        zeroesP = [i for i, x in enumerate(sb1) if x == 0]
        zeroesN = len(zeroesP)

        # putting N-1-nZeroes degrees into N vertices
        for comb2 in itertools.combinations(list(range(2*n - zeroesN - 2)), n-1):
            sb2 = combination_to_sequence(comb2, n - zeroesN - 1)
            for i in zeroesP:
                sb2[i] += 1
            yield tuple(sb1), tuple(sb2)


def directed_rooted_tree_DS(n, bounds=()):
    sb2 = [0 if 0==j else 1 for j in range(n)]
    if bounds:
        for comb in combinations_w_bounds(2*n - 3, n-1, bounds):
            sb = combination_to_sequence(comb, n-2)
            sb1 = [sb[j]+1 if 0 == j else sb[j] for j in range(n)]
            yield tuple(sb1), tuple(sb2)
    else:
        for comb in itertools.combinations(list(range(2*n - 3)), n-1):
            sb = combination_to_sequence(comb, n-2)
            sb1 = [sb[j]+1 if 0 == j else sb[j] for j in range(n)]
            yield tuple(sb1), tuple(sb2)


def combinations_w_bounds(n, r, bounds):
    """https://docs.python.org/3/library/itertools.html#itertools.combinations"""
    indices = [-1] + list(range(r)) + [n]
    # print indices
    for i in reversed(list(range(1, r+1))):
        if indices[i+1] - indices[i] > (bounds[i] + 1):
            indices[i] = indices[i+1] - (bounds[i] + 1)
    # print indices
    yield tuple(indices[1:-1])
    while True:
        # print
        for i in reversed(list(range(1, r+1))):
            if indices[i]+1 - indices[i-1] <= (bounds[i-1] + 1) and indices[i]+1 != indices[i+1]:
                indices[i] += 1
                break
        else:
            return
        # print indices
        for j in range(i+1, r+1):
            indices[j] = indices[j-1] + 1
        # print indices
        for i in reversed(list(range(i+1, r+1))):
            if indices[i+1] - indices[i] > (bounds[i] + 1):
                indices[i] = indices[i+1] - (bounds[i] + 1)

        yield tuple(indices[1:-1])


def combination_to_sequence(array, r):
    M = len(array)
    result = [0 for i in range(M+1)]
    result[0] = array[0]
    for i in range(1, M):
        result[i] = array[i]-array[i-1]-1
    result[M] = r + M - array[M-1] - 1
    return result


if __name__ == '__main__':
    pass
