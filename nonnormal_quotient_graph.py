import sage.all
from nonnormal_quotient import GroupQuotient
from sage.graphs.graph import Graph
from sage.groups.perm_gps.permgroup_named import AlternatingGroup

grp = AlternatingGroup(5)
print("We examine a nonnormal quotient of {}.".format(grp))
print("")

print("In order to choose which quotient to take, we check the order of each subgroup.")
lis = grp.subgroups()
for i in range(len(lis)):
    print([i,len(lis[i])])
print("")

subg = 51
print("In this case, we choose subgroup {}.".format(subg))
if grp.subgroups()[subg].is_normal() is True:
    print("This subgroup is normal in {}.".format(grp))
else:
    print("This subgroup is not normal in {}.".format(grp))
print("")

chir = 'left'
print("We now construct the {} quotient group.".format(chir))
a = GroupQuotient(grp, subg, chir)
print("Below we list the elements of the chosen subgroup, which is {}.".format(a.subgroup))
print(list(a.subgroup))
print("")

print("We create a matrix encoding which blocks indirectly related to the identity block are related to each other.")
mat = a.identity_block_matrix()
for i in mat:
    print(i)
print("")

print("We list the elements of each block with that block's position in the entries of the above matrix.")
tup = tuple(a.identity_related_blocks())
for i in range(len(tup)):
    print([i, a.blocks[tup[i]]])
print("")

print("We produce the graph which has our matrix as its adjacency matrix.")
a.generate_identity_block_relation_graph("relation_graph.svg", layout_choice='spring')
print("")

print("Finally, we examine the associated graphs which identify two vertices in the original graph")
print("if and only if those vertices are minimal distance d apart.")
g = Graph(mat)
for i in range(1,6):
    g.distance_graph(i).plot().save("relation_graph_distance_{}.svg".format(i))
print("")

print("The generated graphs should be saved as svg vector graphics files in the same folder that this program is stored.")
print("By changing the variables 'grp', 'subg', and 'chir' one can see this tutorial run with other examples.")