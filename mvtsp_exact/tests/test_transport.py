import random
from pprint import pprint

from mvtsp_exact.tools.transport import transport

n = 5
random.seed(0)
cost = [[random.random() for i in range(n)] for j in range(n)]
supply = [1, 2, 3, 4, 5]
demand = [5, 4, 3, 2, 1]
pprint(transport(cost, supply, demand))