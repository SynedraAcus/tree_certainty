#! /usr/bin/env python3

import matplotlib.pyplot as plt
import sys

from argparse import ArgumentParser
from dendropy import Tree, TreeList
from dendropy.dataio.newickreader import NewickReader
from methods import *
from sklearn.manifold import MDS
from warnings import catch_warnings, simplefilter

if __name__ == '__main__':
    parser = ArgumentParser('Analysis of consistency in a phylogenetic tree set')
    parser.add_argument('-t', type=str, nargs='*',
                        help='Tree files (multi-Newick). Each file should contain all trees from one analysis')
    parser.add_argument('--unrooted', action='store_true',
                        help='Ignore rooting')
    parser.add_argument('--nonmetric', action='store_true',
                        help='Use nonmetric MDS instead of metric.')
    args = parser.parse_args()
    normalized_shareds = []
    for treefile in args.t:
        print(f'Processing {treefile}...', end='', file=sys.stderr)
        # TODO: change to iterable in case there are too many trees to fit in memory
        try:
            trees = TreeList.get(file=open(treefile), schema='newick')
            if args.unrooted:
                for tree in trees:
                    tree.is_rooted = False
            else:
                # Both values are set explicitly because dendropy reader
                # sets the rootedness in a convoluted way depending on
                # the file, and we want to choose between rooted and unrooted
                # treatment for the entire dataset
                for tree in trees:
                    tree.is_rooted = True
        except NewickReader.NewickReaderDuplicateTaxonError:
            print(' Duplicate taxa names, skipping this file', file=sys.stderr)
            continue
        try:
            vectors = produce_vectors(trees)
        except ValueError as excep:
            # Raised by produce_vectors if the leaf sets differ between trees
            print('\n' + str(excep), file=sys.stderr)
            print('Skipping this file', file=sys.stderr)
            continue
        shared = stable_pair_count(vectors)
        shared_norm = shared / len(vectors[0])
        print(f' {len(vectors)} trees processed, {shared} ({shared_norm * 100:.2f} %) distances not changing',
              file=sys.stderr)
        normalized_shareds.append(shared_norm)
        matrix = vector_distance_matrix(vectors, normalize=True)
        # Produce an image
        with catch_warnings():
            # Disable FutureWarning spam from mds
            simplefilter('ignore')
            mds = MDS(n_components=2, dissimilarity='precomputed', metric=not args.nonmetric)
            transformed = mds.fit_transform(vector_distance_matrix(vectors, normalize=False))
        plt.figure()
        plt.scatter(transformed[:, 0], transformed[:, 1])
        set_name = '.'.join(treefile.split('.')[:-1])
        plt.title(f'Dataset {set_name}, {len(vectors)} trees, {shared} ({shared_norm*100:.2f} %) shared distances')
        plt.savefig(f'{set_name}.png')
        plt.close()
    # Drawing a histogram
    plt.figure()
    plt.hist(normalized_shareds)
    plt.title('Distribution of normalized shared distances')
    plt.savefig('Distances.png')