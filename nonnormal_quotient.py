"""
Nonnormal quotient group tools

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from sage.matrix.constructor import matrix
from sage.graphs.graph import Graph

class GroupQuoitent:
    """
    A quotient object constructed from a group and one of its (possibly nonnormal) subgroups.
    Either left- or right-handed cosets can be used in the construction.

    Attributes:
        group (sage.groups.perm_gps): The initial group used in the construction.
        elems (frozenset): The elements of the initial group.
        subgroup (sage.groups.perm_gps): The chosen subgroup of the initial group.
        cosets (list): The cosets of the chosen subgroup in the initial group.
        comps (list): The components of the relation graph for the relation encoded in `matrix`.
        verbose_comps (list): The same as `comps` but the elements of each entry are group elements rather than indices.
    """

    def __init__(self, group, subgroup, chirality):
        """
        Args:
            group (sage.groups.perm_gps): The initial group used in the construction.
            subgroup (int): The index of a subgroup of that group in the list of subgroups given by ``group.subgroups()``.
            chirality (str): Either 'left' or 'right' for the left or right cosets of that subgroup, respectively.
        """

        self.group = group
        self.elems = tuple(self.group)
        self.subgroup = group.subgroups()[subgroup]
        self.cosets = self.group.cosets(self.subgroup, chirality)
        self.comps = Graph(self.matrix()).connected_components()
        self.verbose_comps = [[self.elems[i] for i in comp] for comp in self.comps]

    def blocks(self):
        """
        Generate blocks in the given group induced by the cosets of the given subgroup.

        Yields:
            frozenset: The elements in the next block.
        """

        for coset_1 in self.cosets:
            for coset_2 in self.cosets:
                block = frozenset(x*y for x in coset_1 for y in coset_2)
                yield block

    def blocks_gens(self):
        """
        Generate a representative pair for each block in the given group induced by the cosets of the given subgroup.

        Yields:
            tuple of sage.groups.perm_gps.permgroup_element.PermutationGroupElement: The representative pair.
        """

        for coset_1 in self.cosets:
            for coset_2 in self.cosets:
                yield (coset_1[0],coset_2[0])

    def block_matrix(self):
        """
        Create a relation matrix for the relation on the blocks induced by the given subgroup.

        Returns:
            sage.matrix.matrix_integer_sparse.Matrix_integer_sparse: A matrix encoding the relation.
        """

        blocks = tuple(self.blocks())
        k = len(blocks)
        m = matrix(k,sparse=True)
        for i in range(k):
            for j in range(i,k):
                if len(set(blocks[i]).intersection(set(blocks[j]))) > 0:
                    m[i,j]=1
                    m[j,i]=1
        return m

    def identity_related_blocks(self):
        """
        Generate the indices assigned to those blocks which are related to the identity block, even indirectly.

        Yields:
            The index of another new block related to the identity block.
        """

        m = self.block_matrix()
        related_blocks = []
        new_blocks = [0]
        while new_blocks != []:
            for i in new_blocks:
                for j in range(len(m)):
                    if m[i,j] == 1 and (j not in related_blocks and j not in new_blocks):
                        new_blocks.append(j)
                        yield j
                related_blocks.append(i)
                new_blocks.remove(i)

    def identity_block_matrix(self):
        """
        Create the submatrix of the block_matrix for those blocks related to the identity, even indirectly.

        Returns:
            sage.matrix.matrix_integer_sparse.Matrix_integer_sparse: The appropriate submatrix of the block_matrix.
        """

        m = self.block_matrix()
        b = tuple(self.identity_related_blocks())
        a = matrix(len(b), sparse=True)
        for i in b:
            for j in b:
                a[b.index(i),b.index(j)] = m[i,j]
        return a

    def generate_identity_block_relation_graph(self, name, layout_choice='spring'):
        """
        Generate the relation graph corresponding to the block relation matrix for those blocks related to the identity block,
        even indirectly.

        ## Fix this
        """

        # the igraph stuff has spotty support right now, couln't get it installed
        # look into GraphData
        # blocks = self.identity_related_blocks
        # labels = [str(self.blocks_gens()[i][0])+","+str(self.blocks_gens()[i][1]) for i in blocks]
        g = Graph(self.identity_block_matrix)
        # g.igraph_graph({'name': labels})
        g.remove_loops()
        g.plot(graph_border=True, layout=layout_choice).save(filename=name)

    def matrix(self):
        """
        Create the relation matrix for the relation on the group elements induced by the given subgroup.

        Returns:
            sage.matrix.matrix_integer_sparse.Matrix_integer_sparse: The relation matrix for the group elements.
        """

        m = matrix(self.group.order(),sparse=True)
        for coset_1 in self.cosets:
            for coset_2 in self.cosets:
                lis = []
                for a in coset_1:
                    for b in coset_2:
                        lis.append(a*b)
                for a in lis:
                    for b in lis:
                        m[self.elems.index(a),self.elems.index(b)]=1
        return m

    def generate_relation_graph(self, name, layout_choice='spring'):
        """
        Generate the relation graph corresponding to the relation matrix on the elements of the group.

        Args:
            name (str): The name of the output file.
            layout_choice (str): The layout format for the graph.
        """

        g = Graph(self.matrix)
        g.remove_loops()
        g.plot(graph_border=True, layout=layout_choice).save(filename=name)

    def group_multiplication_table(self):
        """
        Print a copy of the multiplication table for the given group, up to isomorphism.
        """

        elems = []
        for coset in self.cosets:
            for i in range(len(self.cosets[0])):
                elems.append(coset[i])
        for x in elems:
            row = []
            for y in elems:
                row.append(elems.index(x*y))
            print(row)

    def multiplication_table(self):
        """
        Print the multiplication table for the quotient, up to isomorphism.
        """

        s = [set(comp) for comp in self.verbose_comps]
        for a in s:
            row = []
            for b in s:
                lis = []
                for alpha in a:
                    for beta in b:
                        lis.append(alpha*beta)
                row.append(s.index(set(lis)))
            print(row)

def subgroup_data(group):
    """
    Produce a collection of dictionaries containing data about the subgroups of the given group.

    Args:
        group (sage.groups.perm_gps): The group in question.

    Returns:
        lis of dict: Dictionaries of information about each subgroup of `group`.

    Keys:
        subgroup number - The place the subgroup has in the ordering SageMath assigns to the subgroups of the given group
        subgroup generators - A list of generators for the subgroup
        subgroup normal - True if the subgroup is normal in the given group, Falso otherwise
        quotient group chirality - Left if left cosets were used in the quotient, right otherwise
        quotient group order - The number of elements in the resulting quotient
    """

    lis = []
    for i in range(len(group.subgroups())):
        subgrp = group.subgroups()[i]
        for chirality in ['left', 'right']:
            q = GroupQuoitent(group, i, chirality)
            dic = {}
            dic['subgroup number'] = i
            dic['subgroup generators'] = subgrp.gens()
            dic['subgroup normal'] = subgrp.is_normal()
            dic['quotient group chirality'] = chirality
            dic['quotient group order'] = len(q.comps)
            lis.append(dic)
    return lis

def generate_subgroup_table(group, file_name):
    """
    Write the LaTeX source for a table detailing the information in the subgroup_data dictionary.

    Args:
        group (sage.groups.perm_gps): The group in question.
    """

    lis = subgroup_data(group)
    print("\\begin{tabular}{r | *{4}{c|}}")
    print("subgroup \# & generators & normal? & quotient chirality & quotient order \\\\ \\hline")
    for dic in lis:
        print("%d & %s & %s & %s & %d \\\\ \\cline{2-5}"%(dic['subgroup number'],dic['subgroup generators'],dic['subgroup normal'],dic['quotient group chirality'],dic['quotient group order']))
    print("\end{tabular}")

def normal_closure(group, subgrp_num):
    """
    Return a list of the elements in the normal closure of the given subgroup.
    """
    subgrp = set(list(group.subgroups()[subgrp_num]))
    closure = {}
    for i in range(len(list(group.normal_subgroups()))):
        if subgrp.intersection(set(list(group.normal_subgroups()[i]))) == subgrp:
            if closure == {}:
                closure = set(list(group.normal_subgroups()[i]))
            else:
                closure = closure.intersection(set(list(group.normal_subgroups()[i])))
    return closure

def normal_closure_test(group):
    """
    Check whether the set of elements related, even indirectly, to the identity in the given group under the relation induced by
    taking a quotient is always the normal closure of the subgroup for which the quotient is being taken.  If not, return the 
    elements of the subgroup for which this condition fails to hold.

    This should always return "Test passed".

    Args:
        group: A permutation group.
    """

    for subgrp_num in range(len(group.subgroups())):
        subgroup = group.subgroups()[subgrp_num]
        if subgroup.is_normal() == False:
            if set(GroupQuoitent(group, subgrp_num, 'left').verbose_comps[0]) != normal_closure(group, subgrp_num):
                return(list(group.subgroups()[subgrp_num]))
            if set(GroupQuoitent(group, subgrp_num, 'right').verbose_comps[0]) != normal_closure(group, subgrp_num):
                return(list(group.subgroups()[subgrp_num]))
    print("Test passed.")