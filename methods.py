"""
Methods for measuring the diversity of a set of phylogenetic trees
"""
from dendropy import Tree
from matrix_fix import CorrectMatrix


def distance_vector(tree):
    """
    Produce a distance vector for a given tree
    :param tree:
    :return:
    """
    if not isinstance(tree, Tree):
        raise TypeError(f'{type(tree)} used instead of dendropy.Tree in distance_vector')
    matrix = CorrectMatrix.from_tree(tree)
    # Convert distance matrix into label-sorted vector
    labels = sorted(tree.taxon_namespace.labels())
    distances = []
    for taxon_index, taxon in enumerate(labels):
        for other_taxon in labels[taxon_index + 1:]:
            distances.append(matrix.path_edge_count(tree.taxon_namespace.get_taxon(taxon),
                                                    tree.taxon_namespace.get_taxon(other_taxon)) - 1)
    return distances


def produce_vectors(trees):
    """
    Produce vectors for all trees in the iterable.

    This function modifies trees in the process, so use a copy of the tree set
    if you expect to need it later.
    :param trees:
    :return:
    """
    orig_labels = trees[0].taxon_namespace.labels()
    count = 0
    vectors = []
    for tree in trees:
        count += 1
        labels = tree.taxon_namespace.labels()
        if labels != orig_labels:
            raise ValueError(f'Label set in tree #{count} is different from previous trees')
        vectors.append(distance_vector(tree))
    return vectors


def shared_pairs(vector1: list[int], vector2: list[int],
                 normalize: bool = False):
    """
    Take two distance vectors, return the count of taxa pairs that have the same
    patristic distances in both.
    :param vector1: distance vector
    :param vector2: distance vector
    :param normalize: Whether the count should be normalized (divided by vector length)
    . Defaults to False, ie returning raw counts
    :return:
    """
    if len(vector1) != len(vector2):
        raise ValueError('Vector lengths do not match')
    count = sum(map(lambda x: 1 if x[0] == x[1] else 0, zip(vector1, vector2)))
    if normalize:
        # TODO: maybe other normalization procedures besides simple division
        return count / len(vector1)
    return count


def tree_similarity_matrix(trees: list[Tree], **kwargs):
    """
    For an list of trees, return a similarity matrix
    :param trees:
    :return:
    """
    # Does not check tree validity, expects it to be checked by caller
    # Or raised by shared_pairs
    vectors = [distance_vector(x) for x in trees]
    return vector_similarity_matrix(vectors, **kwargs)


def vector_similarity_matrix(vectors: list[list[int]], **kwargs):
    """
    Build a similarity matrix for a list of distance vectors
    :param vectors:
    :param kwargs:
    :return:
    """
    matrix = [[0 for _ in vectors] for _1 in vectors]
    for index, vector in enumerate(vectors):
        matrix[index][index] = 0
        for other_index, other_vector in enumerate(vectors[index:]):
            dist = shared_pairs(vector, other_vector, **kwargs)
            # Does an unnecessary assignment when other_index = 0, but who cares
            matrix[index][index + other_index] = dist
            matrix[index + other_index][index] = dist
    return matrix


def vector_distance_matrix(vectors, normalize=False):
    """
    Build a distance matrix for a list of distance vectors
    :param vectors:
    :param normalize:
    :return:
    """
    matrix = vector_similarity_matrix(vectors, normalize = normalize)
    r = []
    if normalize:
        # If normalized, distance = 1 - similarity
        for row in matrix:
            r.append([1 - x for x in row])
    else:
        for row in matrix:
            max_value = len(vectors[0])
            r.append([max_value - x for x in row])
    return r


def tree_distance_matrix(trees, **kwargs):
    """
    Build a distance matrix for a list of trees
    :param trees:
    :param kwargs:
    :return:
    """
    vectors = [distance_vector(x) for x in trees]
    return vector_distance_matrix(vectors, **kwargs)


def stable_pairs(vectors):
    """
    Take an iterable of distance vectors, return the set of taxa pairs that
    have the same patristic distance in all vectors.
    This function does not calculate a complete distance matrix and is
    intended to identify the stable taxa pairs relatively quickly in a large
    set of vectors. If you need pairwise comparisons, use dedicated methods
    (`tree_distance_matrix` or `vector_distance_matrix`) instead.
    :param vectors:
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


def stable_pair_count(vectors, normalize=False):
    """
    A number of stable pairs, either normalized or not
    :param vectors:
    :param normalize:
    :return:
    """
    is_stable = stable_pairs(vectors)
    count = is_stable.count(True)
    if normalize:
        return count / len(vectors[0])
    else:
        return count