import sage.all
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from nonnormal_quotient import normal_closure, normal_closure_test

n = 5
print("Let g be the symmetric group on {} elements.".format(n))
g = SymmetricGroup(n)
print(g)
print("")
print("Pick a subgroup of g.")
print(g.subgroups()[1])
print("")
print("List the elements of the normal closure of that subgroup in g.")
print(list(normal_closure(g,1)))
normal_closure_test(g)