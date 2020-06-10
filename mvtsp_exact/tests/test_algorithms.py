import time
from pprint import pprint
import random

from mvtsp_exact.dc_mv import dc_MV
from mvtsp_exact.dp_mv import dp_MV
from mvtsp_exact.enum_ds import enum_DS
from mvtsp_exact.enum_mv import enum_MV

n = 5
random.seed(0)
cost = [[random.random() for i in range(n)] for j in range(n)]
r = [random.randint(1, 10) for i in range(n)]
pprint(cost, width=160)
print(r)
print()

print("enum-MV")
t = time.time()
pprint(enum_MV(cost, r), width=140)
print("{:.4f} s\n".format(time.time()-t))

print()

print("enum-DS")
t = time.time()
pprint(enum_DS(cost, r), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("enum-DS, rooted")
t = time.time()
pprint(enum_DS(cost, r, rootedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("enum-DS, rooted, bounded")
t = time.time()
pprint(enum_DS(cost, r, rootedQ=True, boundedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print()

print("dp-MV")
t = time.time()
pprint(dp_MV(cost, r), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dp-MV, rooted")
t = time.time()
pprint(dp_MV(cost, r, rootedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dp-MV, rooted, bounded")
t = time.time()
pprint(dp_MV(cost, r, rootedQ=True, boundedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print()

print("dc-MV")
t = time.time()
pprint(dc_MV(cost, r), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dc-MV, rooted")
t = time.time()
pprint(dc_MV(cost, r, rootedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dc-MV, rooted, bounded")
t = time.time()
pprint(dc_MV(cost, r, rootedQ=True, boundedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dc-MV, perfectly balanced, rooted")
t = time.time()
pprint(dc_MV(cost, r, perfectlyBalancedQ=True, rootedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))

print("dc-MV, perfectly balanced, rooted, bounded")
t = time.time()
pprint(dc_MV(cost, r, perfectlyBalancedQ=True, rootedQ=True, boundedQ=True), width=140)
print("{:.4f} s\n".format(time.time()-t))
