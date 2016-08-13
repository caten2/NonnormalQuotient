"""
Nonnormal quotient closure conjecture example

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from nonnormal_quotient import *

n = 4
print('Let `g` be the symmetric group on {} elements.'.format(n))
g = SymmetricGroup(n)
print(g)
print('')

print('Pick a subgroup of `g`.')
print(g.subgroups()[1])
print('')

print('List the elements of the normal closure of that subgroup in `g`.')
print(list(find_normal_closure(g,1)))
print('')

print('Test whether the equivalence class of the identity block is the normal closure of the given subgroup for any subgroup of `g`.')
perform_normal_closure_test(g)