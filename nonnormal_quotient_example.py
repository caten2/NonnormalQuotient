"""
Nonnormal quotient example

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from nonnormal_quotient import GroupQuotient, tex_subgroup_table

print('Let `g` be the symmetric group on 4 elements.')
g = SymmetricGroup(4)
print('Generate LaTeX source for a table showing information about the subgroups of g.')
tex_subgroup_table(g)
print('')

print('Take the left-handed quotient by subgroup 1, which is not normal in g.')
k = GroupQuotient(g,1,'left')
print('Display the sets of related elements.')
print(k.verbose_components)
print('')

print('Save the relation graph as a png image file called \'relation_graph.png\'.')
k.draw_relation_graph('relation_graph.png')
print('')

print('Display the multiplication table for this quotient.')
k.display_multiplication_table()
print('')

print('Generate LaTeX markup for the detailed multiplication table for this quotient.')
k.tex_detailed_multiplication_table()