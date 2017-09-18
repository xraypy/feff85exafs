import numpy as np
from feffpath import ScatteringPath

path = ScatteringPath(phase_file='phase.pad')
assert(path is not None)

path.set_absorber( x=0.01,   y=0.1,   z=0.01)
path.add_scatterer(x=1.8058, y=0.005, z=1.8063, ipot=1)
path.degen = 12


path.calculate_xafs()

assert(path.nleg == 2)
assert(path.iz[0] == 29)
assert(path.rat[0, 0] > 0.00)
assert(path.rat[0, 0] < 0.20)

assert(path.rat[0, 1] > 1.70)
assert(path.rat[0, 1] < 1.90)

assert(path.rnorman > 2.50)
assert(path.rnorman < 2.80)

assert(path.rs   > 1.80)
assert(path.rs   < 2.10)
assert(path.kf   > 0.80)
assert(path.kf   < 1.10)

assert(path.exch_label.startswith('H-L'))


npts = 1 + max(np.where(path.kfeff > 0)[0])
assert(npts > 25)

assert(path.kfeff[5] > 0.490)
assert(path.kfeff[5] < 0.510)

assert(path.rep[5] > 1.880)
assert(path.rep[5] < 1.900)

assert(path.lam[5] > 16.4)
assert(path.lam[5] < 16.6)

assert(path.real_phc[5] > 3.92)
assert(path.real_phc[5] < 3.98)
