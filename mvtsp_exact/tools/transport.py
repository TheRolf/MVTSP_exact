import networkx as nx


def transport(cost, supply, demand):
    N = len(cost)
    G = nx.DiGraph()
    G.add_node('s', demand=-sum(supply))
    G.add_node('t', demand=sum(demand))
    # G.add_edge('t', 's', weight=0, capacity=sum(r))
    for i in range(N):
        G.add_edge('s', 's_' + str(i), weight=0, capacity=supply[i])
        G.add_edge('t_' + str(i), 't', weight=0, capacity=demand[i])
        for j in range(N):
            if cost[i][j] < float("inf"):
                G.add_edge('s_' + str(i), 't_' + str(j), weight=cost[i][j], capacity=float('inf'))
    flowCost, flowDict = nx.capacity_scaling(G)
    flowXprime = {}
    for u in flowDict:
        for v in flowDict[u]:
            if "_" in u and "_" in v and flowDict[u][v] > 0:
                i, j = int(u[2:]), int(v[2:])
                flowXprime[(i, j)] = flowDict[u][v]
    return flowCost, flowXprime


if __name__ == '__main__':
    pass
