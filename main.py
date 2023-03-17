import matplotlib.pyplot as plt

from argparse import ArgumentParser
from dendropy import Tree, TreeList
from methods import *
from sklearn.manifold import MDS

if __name__ == '__main__':
    parser = ArgumentParser('Analysis of consistency in a phylogenetic tree set')
    parser.add_argument('-t', type=str, help='Tree file (multi-Newick)')
    parser.add_argument('--unrooted', action='store_true',
                        help='Ignore rooting')
    args = parser.parse_args()
    # TODO: work with unrooted trees correctly
    # Dendropy implementation counts the root anyway, so the distance between
    # leaves on different halves of the tree is overestimated.
    if args.unrooted:
        raise NotImplementedError('Unrooted tree analysis is not ready yet')
    trees = TreeList.get(file=open(args.t), schema='newick')
    # TODO: change to iterable in case there are too many trees to fit in memory
    # Validating and modifying trees: check that all share the same leaf labels, set all branch lengths to 1.0
    orig_labels = trees[0].taxon_namespace.labels()
    count = 0
    vectors = []
    for tree in trees:
        count += 1
        labels = tree.taxon_namespace.labels()
        tree.is_rooted = True
        if labels != orig_labels:
            raise ValueError(f'Label set in tree #{count} is different from previous trees')
        vectors.append(distance_vector(tree))
    shared = stable_pair_count(vectors)
    shared_norm = shared / len(vectors[0])
    matrix = vector_distance_matrix(vectors, normalize=True)
    # Produce an image
    # TODO: multiple tree files, multiple images
    mds = MDS(n_components=2, dissimilarity='precomputed')
    transformed = mds.fit_transform(vector_distance_matrix(vectors, normalize=False))
    plt.scatter(transformed[:, 0], transformed[:, 1])
    set_name = '.'.join(args.t.split('.')[:-1])
    plt.title(f'File {set_name}, {shared} ({shared_norm*100:.2f} %) shared distances')
    plt.savefig(f'{set_name}.png')
