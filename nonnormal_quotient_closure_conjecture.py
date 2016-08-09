import sage.all
from sage.groups.perm_gps.permgroup_named import SymmetricGroup,\
    AlternatingGroup, DihedralGroup
from nonnormal_quotient import *
from sage.graphs.graph import Graph

# n = 4
# print('Let g be the symmetric group on {} elements.'.format(n))
# g = SymmetricGroup(n)
# print(g)
# print('')
# print('Pick a subgroup of g.')
# print(g.subgroups()[1])
# print('')
# print('List the elements of the normal closure of that subgroup in g.')
# print(list(normal_closure(g,1)))
# print('')
# normal_closure_test(g)

# group = SymmetricGroup(4)
# # identity_block_transitivity_test(group)
# transitivity_test(group)

quot = GroupQuotient(AlternatingGroup(5),1,'left')
print(relation_matrix_transitivity_test(quot.block_matrix()))
# g = Graph(quot.matrix())
# g = g.connected_components_subgraphs()[0]
# g.plot(size=[100,100]).save('asdf.svg')