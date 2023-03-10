"""
Methods for measuring the diversity of a set of phylogenetic trees
"""
from dendropy import Tree


def distance_vector(tree):
    """
    Produce a distance vector for a given tree
    :param tree:
    :return:
    """
    if not isinstance(tree, Tree):
        raise TypeError(f'{type(tree)} used instead of dendropy.Tree in distance_vector')
    matrix = tree.phylogenetic_distance_matrix()
    # Convert distance matrix into label-sorted vector
    labels = sorted(tree.taxon_namespace.labels())
    distances = []
    for taxon_index, taxon in enumerate(labels):
        for other_taxon in labels[taxon_index + 1:]:
            distances.append(matrix.path_edge_count(tree.taxon_namespace.get_taxon(taxon),
                                                    tree.taxon_namespace.get_taxon(other_taxon)) - 1)
            print(taxon, other_taxon, distances[-1])
    return distances


def stable_pairs(vectors):
    """
    Take an iterable of matrices, return the set of taxa pairs that have the
    same patristic distance in all matrices
    :param matrices:
    :return:
    """
    is_stable = None
    for vector in vectors:
        if not is_stable:
            # On the first tree, set the vector collection and results list
            is_stable = [True for _ in vector]
            orig_vector = vector
            continue
        if len(vector) != len(is_stable):
            raise ValueError('Incorrect vector length')
        for index in range(len(vector)):
            if is_stable[index] and vector[index] != orig_vector[index]:
                is_stable[index] = False
    return is_stable