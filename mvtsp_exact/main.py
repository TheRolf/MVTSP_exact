import random
from pprint import pprint

from mvtsp_exact.dc_mv import dc_MV
from mvtsp_exact.dp_mv import dp_MV
from mvtsp_exact.enum_ds import enum_DS
from mvtsp_exact.enum_mv import enum_MV

# random instance, setting r=1 for some cities for more interesting cases
n = 5
r = [1 if random.random() < 0.5 else random.randint(1, 100) for i in range(n)]
c = [[random.random() for i in range(n)] for j in range(n)]

print("r =", r)
print("c = ", end="")
pprint(c, width=300)

# dp_MV (fastest)
print("\ndp_MV")
cost, solution = dp_MV(c, r)
print(cost, solution)

# dc_MV
print("\ndc_MV")
cost, solution = dc_MV(c, r)
print(cost, solution)

# enum_DS
print("\nenum_DS")
cost, solution = enum_DS(c, r)
print(cost, solution)

# enum_MV (takes a looong time for n>5)
print("\nenum_MV")
cost, solution = enum_MV(c, r)
print(cost, solution)
