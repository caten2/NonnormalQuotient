"""
Nonnormal quotient group tools

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from sage.matrix.constructor import matrix
from sage.graphs.graph import Graph

class GroupQuotient:
    """
    A quotient object constructed from a group and one of its (possibly nonnormal) subgroups.
    Either left- or right-handed cosets can be used in the construction.

    Attributes:
        group (sage.groups.perm_gps): The initial group used in the construction.
        elements (frozenset): The elements of the initial group.
        subgroup (sage.groups.perm_gps): The chosen subgroup of the initial group.
        chirality (str): Either 'left' or 'right' for the left or right cosets of that subgroup, respectively.
        cosets (tuple): The cosets of the chosen subgroup in the initial group.
        blocks (tuple): The blocks of the chosen subgroup in the initial group.
        components (tuple): The components of the relation graph for the relation encoded in ``self.create_relation_matrix()``.
        verbose_components (tuple): The same as `self.components` but the elements of each entry are group elements rather than indices.
    """

    def __init__(self, group, subgroup, chirality):
        """
        Args:
            group (sage.groups.perm_gps): The initial group used in the construction.
            subgroup (int): The index of a subgroup of that group in the list of subgroups given by ``group.subgroups()``.
            chirality (str): Either 'left' or 'right' for the left or right cosets of that subgroup, respectively.
        """

        self.group = group
        self.elements = tuple(self.group)
        self.subgroup = group.subgroups()[subgroup]
        self.chirality = chirality
        self.cosets = tuple(self.group.cosets(self.subgroup,chirality))
        self.blocks = tuple(self.generate_blocks())
        self.components = tuple(Graph(self.create_relation_matrix()).connected_components())
        self.verbose_components = tuple([self.elements[i] for i in comp] for comp in self.components)

    def generate_blocks(self):
        """
        Generate blocks in the given group induced by the cosets of the given subgroup.

        Yields:
            frozenset: The elements in the next block.
        """

        for coset_1 in self.cosets:
            for coset_2 in self.cosets:
                block = frozenset(x*y for x in coset_1 for y in coset_2)
                yield block

    def generate_representative_pairs(self):
        """
        Generate a representative pair for each block in the given group induced by the cosets of the given subgroup.

        Yields:
            tuple of sage.groups.perm_gps.permgroup_element.PermutationGroupElement: The representative pair.
        """

        for coset_1 in self.cosets:
            for coset_2 in self.cosets:
                yield (coset_1[0],coset_2[0])

    def create_block_matrix(self):
        """
        Create a relation create_relation_matrix for the relation on the blocks induced by the given subgroup. This relation is called '$\rho$' in the associated paper.

        Returns:
            sage.create_relation_matrix.matrix_integer_sparse.Matrix_integer_sparse: A create_relation_matrix encoding the relation.
        """

        k = len(self.blocks)
        m = matrix(k,sparse=True)
        for i in range(k):
            for j in range(i,k):
                if set(self.blocks[i]).intersection(set(self.blocks[j])) != set():
                    m[i,j]=1
                    m[j,i]=1
        return m

    def generate_identity_related_blocks(self):
        """
        Generate the indices assigned to those blocks which are related to the identity block, even indirectly.

        Yields:
            int: The index of another new block related to the identity block.
        """

        m = self.create_block_matrix()
        related_blocks = []
        new_blocks = [0]
        yield 0
        while new_blocks != []:
            for i in new_blocks:
                for j in range(m.dimensions()[0]):
                    if m[i,j] == 1 and (j not in related_blocks and j not in new_blocks):
                        new_blocks.append(j)
                        yield j
                related_blocks.append(i)
                new_blocks.remove(i)

    def create_identity_block_matrix(self):
        """
        Create the submatrix of the create_block_matrix for those blocks related to the identity, even indirectly.

        Returns:
            sage.create_relation_matrix.matrix_integer_sparse.Matrix_integer_sparse: The appropriate submatrix of the create_block_matrix.
        """

        m = self.create_block_matrix()
        b = tuple(self.generate_identity_related_blocks())
        a = matrix(len(b), sparse=True)
        for i in b:
            for j in b:
                a[b.index(i),b.index(j)] = m[i,j]
        return a

    def draw_identity_block_relation_graph(self, name='identity_block_relation_graph.png', layout_choice='spring'):
        """
        Draw the relation graph corresponding to the block relation create_relation_matrix for those blocks related to the identity block,
        even indirectly.

        Args:
            name (str): The name of the output file, optionally with extension.
            layout_choice (str): The layout algorthm used to draw the graph.
        """

        g = Graph(self.create_identity_block_matrix())
        g.remove_loops()
        g.plot(graph_border=False,layout=layout_choice).save(filename=name)

    def create_relation_matrix(self):
        """
        Create the relation create_relation_matrix for the relation on the group elements induced by the given subgroup. This relation is called '$\psi$' in the associated paper.

        Returns:
            sage.create_relation_matrix.matrix_integer_sparse.Matrix_integer_sparse: The relation create_relation_matrix for the group elements.
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
                        m[self.elements.index(a),self.elements.index(b)]=1
        return m

    def draw_relation_graph(self, name, layout_choice='spring'):
        """
        Draw the relation graph corresponding to the relation create_relation_matrix on the elements of the group.

        Args:
            name (str): The name of the output file.
            layout_choice (str): The layout format for the graph.
        """

        g = Graph(self.create_relation_matrix())
        g.remove_loops()
        g.plot(graph_border=False, layout=layout_choice).save(filename=name)

    def draw_block_relation_graph(self, name, layout_choice='spring'):
        """
        Draw the relation graph corresponding to the relation create_block_matrix on the elements of the group.

        Args:
            name (str): The name of the output file.
            layout_choice (str): The layout format for the graph.
        """

        g = Graph(self.create_block_matrix())
        g.remove_loops()
        g.plot(graph_border=False, layout=layout_choice).save(filename=name)

    def display_group_multiplication_table(self):
        """
        Print a copy of the multiplication table for the given group, up to isomorphism.
        """

        elems = []
        for coset in self.cosets:
            for i in range(len(self.cosets[0])):
                elems.append(coset[i])
        elems = tuple(elems)
        for x in elems:
            row = []
            for y in elems:
                row.append(elems.index(x*y))
            print(row)

    def display_multiplication_table(self):
        """
        Print the multiplication table for the quotient, up to isomorphism.
        """

        s = tuple(set(comp) for comp in self.verbose_components)
        for a in s:
            row = []
            for b in s:
                lis = []
                for alpha in a:
                    for beta in b:
                        lis.append(alpha*beta)
                row.append(s.index(set(lis)))
            print(row)

    def tex_detailed_multiplication_table(self, file_name='multiplication_table.tex'):
        """
        Write the LaTeX markup for a table showing multiplication in the nonnormal quotient and the corresponding closure.

        Args:
            file_name (str): The name of the file to be generated.
        """

        size = self.group.order()
        sorted_elems = []
        norm_closure = find_normal_closure(self.group,self.group.subgroups().index(self.subgroup),output='group')
        closure_cosets = tuple(self.group.cosets(norm_closure))
        cosets = self.group.cosets(self.subgroup,self.chirality)
        m = len(closure_cosets[0])
        n = len(cosets[0])
        for closure_coset in closure_cosets:
            for coset in cosets:
                intersection = frozenset(closure_coset).intersection(frozenset(coset))
                if intersection != frozenset():
                    sorted_elems.append(tuple(intersection))
        file = open(file_name, 'w')
        file.write('\\documentclass{standalone}\n')
        file.write('\\usepackage{amsmath}\n')
        file.write('\\usepackage{multirow}\n\n')
        file.write('\\providecommand{\\nc}[1]{\\operatorname{Nc}(#1)}\n\n')
        file.write('\\begin{document}\n')
        file.write('\\begin{tabular}{| *{3}{r|} ')
        t = len(cosets)
        file.write(('*{{{}}}{{c}} | '.format(size/t))*t)
        file.write('}} \\cline{{4-{}}}\n'.format(size+3))
        # normal closure cosets
        file.write('\\multicolumn{3}{c|}{} & ')
        s = len(closure_cosets[0])
        for closure_coset in closure_cosets[:-1]:
            rep = ''
            if closure_coset[0] != self.group[0]:
                rep = closure_coset[0]
            file.write('\multicolumn{{{}}}{{c|}}{{${}\\nc{{H}}$}} & '.format(s,rep))
        rep = ''
        if closure_cosets[-1][0] != self.group[0]:
            rep = closure_cosets[-1][0]
        file.write('\multicolumn{{{}}}{{c|}}{{${}\\nc{{H}}$}} '.format(s,rep))
        file.write('\\\\ \\cline{{4-{}}}\n'.format(size+3))
        # nonnormal subgroup cosets
        file.write('\\multicolumn{3}{c|}{} & ')
        t = len(cosets[0])
        for coset in cosets[:-1]:
            rep = ''
            if coset[0] != self.group[0]:
                rep = coset[0]
            file.write('\multicolumn{{{}}}{{c|}}{{${}H$}} & '.format(t,rep))
        rep = ''
        if cosets[-1][0] != self.group[0]:
            rep = cosets[-1][0]
        file.write('\multicolumn{{{}}}{{c|}}{{${}H$}} '.format(t,rep))
        file.write('\\\\ \\cline{{4-{}}}\n'.format(size+3))
        # elements row
        file.write('\\multicolumn{3}{c|}{} & ')
        for i in range(size-1):
            file.write('${}$ & '.format(self.cosets[i/n][i%n]))
        file.write('${}$ '.format(self.cosets[-1][-1]))
        file.write('\\\\ \\hline\n')
        # table values
        for i in range(size):
            if i % m == 0:
                rep = ''
                if closure_cosets[i/m][0] != self.group[0]:
                    rep = closure_cosets[i/m][0]
                file.write('\multirow{{{}}}{{*}}{{${}\\nc{{H}}$}} & '.format(s,rep))
            else:
                file.write('& ')
            if i % n == 0:
                rep = ''
                if cosets[i/n][0] != self.group[0]:
                    rep = cosets[i/n][0]
                file.write('\multirow{{{}}}{{*}}{{${}H$}} & '.format(t,rep))
            else:
                file.write('& ')
            file.write('${}$ & '.format(self.cosets[i/n][i%n]))
            for j in range(size-1):
                file.write('${}$ & '.format(self.cosets[i/n][i%n]*self.cosets[j/n][j%n]))
            file.write('${}$ \\\\'.format(self.cosets[i/n][i%n]*self.cosets[-1][-1]))
            if (i+1) % n == 0:
                if (i+1) % m == 0:
                    file.write(' \\hline\n')
                else:
                    file.write(' \\cline{{2-{}}}\n'.format(size+3))
            else:
                file.write('\n')
        file.write('\end{tabular}\n')
        file.write('\end{document}')
        file.close()

def compute_subgroup_data(group):
    """
    Produce a collection of dictionaries containing data about the subgroups of the given group.

    Args:
        group (sage.groups.perm_gps): The group in question.

    Returns:
        lis of dict: Dictionaries of information about each subgroup of `group`.

    Keys:
        subgroup number - The place the subgroup has in the ordering SageMath assigns to the subgroups of the given group.
        subgroup generators - A list of generators for the subgroup.
        subgroup normal - True if the subgroup is normal in the given group, Falso otherwise.
        quotient group chirality - Left if left cosets were used in the quotient, right otherwise.
        quotient group order - The number of elements in the resulting quotient.
    """

    lis = []
    for i in range(len(group.subgroups())):
        subgrp = group.subgroups()[i]
        for chirality in ['left', 'right']:
            q = GroupQuotient(group, i, chirality)
            dic = {}
            dic['subgroup number'] = i
            dic['subgroup generators'] = subgrp.gens()
            dic['subgroup normal'] = subgrp.is_normal()
            dic['quotient group chirality'] = chirality
            dic['quotient group order'] = len(q.components)
            lis.append(dic)
    return lis

def tex_subgroup_table(group, file_name='subgroup_table.tex'):
    """
    Write the LaTeX source for a table detailing the information in the compute_subgroup_data dictionary.

    Args:
        group (sage.groups.perm_gps): The group in question.
        file_name (str): The name of the file to be generated.
    """

    lis = compute_subgroup_data(group)
    file = open(file_name, 'w')
    file.write('\\documentclass{standalone}\n\n')
    file.write('\\begin{document}\n')
    file.write('\\begin{tabular}{r | *{3}{c|} c}\n')
    file.write('subgroup number & generators & normal & quotient chirality & quotient order \\\\ \\hline\n')
    for dic in lis:
        if dic['subgroup number'] == len(group.subgroups())-1 and dic['quotient group chirality'] == 'right':
            file.write('{} & {} & {} & {} & {} \\\\\n'.format(dic['subgroup number'],', '.join(repr(gen) for gen in dic['subgroup generators']),dic['subgroup normal'],dic['quotient group chirality'],dic['quotient group order']))
        else:
            file.write('{} & {} & {} & {} & {} \\\\ \\hline\n'.format(dic['subgroup number'],', '.join(repr(gen) for gen in dic['subgroup generators']),dic['subgroup normal'],dic['quotient group chirality'],dic['quotient group order']))
    file.write('\end{tabular}\n')
    file.write('\end{document}')
    file.close()

def find_normal_closure(group, subgrp_num, output='set'):
    """
    Return the normal closure of the given subgroup.

    Args:
        group (sage.groups.perm_gps): The group in question.
        subgrou_num (int): The place the subgroup has in the ordering SageMath assigns to the subgroups of the given group.
        output (str): Returns a set if set to 'set' and returns a subgroup of `group` if set to 'group'.

    Returns:
        set: The elements in the normal closure.
        sage.groups.perm_gps: The normal closure as a subgroup of `group`.
    """

    subgrp = frozenset(group.subgroups()[subgrp_num])
    closure = set()
    for i in range(len(list(group.normal_subgroups()))):
        if subgrp.intersection(frozenset(group.normal_subgroups()[i])) == subgrp:
            if closure == set():
                closure = set(list(group.normal_subgroups()[i]))
            else:
                closure = closure.intersection(set(group.normal_subgroups()[i]))
    if output == 'set':
        return closure
    if output == 'group':
        for subg in group.subgroups():
            if set(subg) == closure:
                return subg

def perform_normal_closure_test(group):
    """
    Check whether the set of elements related, even indirectly, to the identity in the given group under the relation induced by
    taking a quotient is always the normal closure of the subgroup for which the quotient is being taken. If not, return the 
    elements of the subgroup for which this condition fails to hold.

    This should always print 'Test passed.'.

    Args:
        group: A permutation group.

    Returns:
        list: The elements of the subgroup for which the test fails.
    """

    for subgrp_num in range(len(group.subgroups())):
        subgroup = group.subgroups()[subgrp_num]
        if subgroup.is_normal() == False:
            if set(GroupQuotient(group, subgrp_num, 'left').verbose_components[0]) != find_normal_closure(group, subgrp_num):
                return(list(group.subgroups()[subgrp_num]))
            if set(GroupQuotient(group, subgrp_num, 'right').verbose_components[0]) != find_normal_closure(group, subgrp_num):
                return(list(group.subgroups()[subgrp_num]))
    print('Test passed.')

def perform_relation_matrix_transitivity_test(mat):
    """
    Check whether a square create_relation_matrix encodes a transitive relation.

    Args:
        mat (create_relation_matrix): The create_relation_matrix in question.

    Returns:
        tuple of int: A triple (i,j,k) that witnesses the failure of `mat` to be transitive.
    """

    size = mat.dimensions()[0]
    for i in range(size):
        for j in range(i,size):
            for k in range(j,size):
                if mat[i][j] == mat[j][k] == 1 and mat[i][k] == 0:
                    return((i,j,k))

def perform_block_transitivity_test(group):
    """
    Test every subgroup of a given group for a transitive block relation.

    Args:
        group (sage.groups.perm_gps): The group in question.
    """

    for subgroup_num in range(len(group.subgroups())):
        quot = GroupQuotient(group,subgroup_num,'left')
        mat = quot.create_block_matrix()
        tup = perform_relation_matrix_transitivity_test(mat)
        if tup is not None:
            i,j,k = tup[0],tup[1],tup[2]
            print(subgroup_num)
            print(quot.subgroup)
            print([quot.blocks[i],quot.blocks[j],quot.blocks[k]])
            return
    print("No counterexample found.")

def perform_identity_block_transitivity_test(group):
    """
    Test whether the restriction of the block relation to those blocks related, even indirectly, to the identity block is transitive.

    Args:
        group (sage.groups.perm_gps): The group in question.
    """

    for subgroup_num in range(len(group.subgroups())):
        quot = GroupQuotient(group,subgroup_num,'left')
        mat = quot.create_identity_block_matrix()
        tup = perform_relation_matrix_transitivity_test(mat)
        if tup is not None:
            i,j,k = tup[0],tup[1],tup[2]
            print(subgroup_num)
            print(quot.subgroup)
            return
    print("No counterexample found.")

def perform_transitivity_test(group):
    """
    For each subgroup of the given group test whether the relation on group elements induced by the corresponding blocks is transitive.

    Args:
        group (sage.groups.perm_gps): The group in question.
    """

    for subgroup_num in range(len(group.subgroups())):
        quot = GroupQuotient(group,subgroup_num,'left')
        mat = quot.create_relation_matrix()
        tup = perform_relation_matrix_transitivity_test(mat)
        if tup is not None:
            i,j,k = tup[0],tup[1],tup[2]
            print(subgroup_num)
            print(quot.subgroup)
            print([quot.elements[i],quot.elements[j],quot.elements[k]])
            return
    print("No counterexample found.")