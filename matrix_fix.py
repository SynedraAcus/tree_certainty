from dendropy import PhylogeneticDistanceMatrix, Tree


class CorrectMatrix(PhylogeneticDistanceMatrix):
    """
    A class that overrides native phylogenetic matrix calculation to correctly
    process distances across the root.
    """
    def compile_from_tree(self, tree):
        """
        Calculate the pairwise distances across the tree.
        """
        self.clear()
        self.taxon_namespace = tree.taxon_namespace
        self._tree_length = 0.0
        self._num_edges = 0
        if self.is_store_path_edges:
            default_pedges = []
        else:
            default_pedges = None
        for node in tree.postorder_node_iter():
            try:
                self._tree_length += node.edge.length
            except TypeError:  # None for edge length
                pass
            self._num_edges += 1
            children = node.child_nodes()
            if len(children) == 0:
                node.desc_paths = {node: (0, 0, default_pedges)}
            else:
                node.desc_paths = {}
                for cidx1, c1 in enumerate(children):
                    for desc1, (desc1_plen, desc1_psteps, desc1_pedges) in c1.desc_paths.items():
                        if c1.edge_length is None:
                            c1_edge_length = 0.0
                        else:
                            c1_edge_length = c1.edge.length
                        if self.is_store_path_edges:
                            pedges = list(desc1_pedges + [c1.edge])
                        else:
                            pedges = default_pedges
                        node.desc_paths[desc1] = (desc1_plen + c1_edge_length, desc1_psteps + 1, pedges)
                        assert desc1.taxon is not None
                        if desc1.taxon not in self._taxon_phylogenetic_distances:
                            self._mapped_taxa.add(desc1.taxon)
                            self._taxon_phylogenetic_distances[desc1.taxon] = {}
                            self._taxon_phylogenetic_distances[desc1.taxon][desc1.taxon] = 0.0
                            self._taxon_phylogenetic_path_steps[desc1.taxon] = {}
                            self._taxon_phylogenetic_path_steps[desc1.taxon][desc1.taxon] = 0
                            if self.is_store_path_edges:
                                self._taxon_phylogenetic_path_edges[desc1.taxon] = {}
                                self._taxon_phylogenetic_path_edges[desc1.taxon][desc1.taxon] = []
                            self._mrca[desc1.taxon] = {desc1.taxon: desc1}
                        for c2 in children[cidx1 + 1:]:
                            for desc2, (desc2_plen, desc2_psteps, desc2_pedges) in c2.desc_paths.items():
                                self._mapped_taxa.add(desc2.taxon)
                                self._mrca[desc1.taxon][desc2.taxon] = c1.parent_node
                                # self._all_distinct_mapped_taxa_pairs.add( tuple([desc1.taxon, desc2.taxon]) )
                                self._all_distinct_mapped_taxa_pairs.add(frozenset([desc1.taxon, desc2.taxon]))
                                if c2.edge_length is None:
                                    c2_edge_length = 0.0
                                else:
                                    c2_edge_length = c2.edge.length
                                pat_dist = node.desc_paths[desc1][0] + desc2_plen + c2_edge_length
                                self._taxon_phylogenetic_distances[desc1.taxon][desc2.taxon] = pat_dist
                                # Correct calculation of step number across the root
                                # Note that it introduces disrepancy between _taxon_phylogenetic_path_edges and
                                # _taxon_phylogenetic_path_steps
                                if not tree.is_rooted and node is tree.seed_node:
                                    path_steps = node.desc_paths[desc1][1] + desc2_psteps
                                else:
                                    path_steps = node.desc_paths[desc1][1] + desc2_psteps + 1
                                self._taxon_phylogenetic_path_steps[desc1.taxon][desc2.taxon] = path_steps
                                if self.is_store_path_edges:
                                    pedges = tuple(node.desc_paths[desc1][2] + [c2.edge] + desc2_pedges[::-1])
                                    self._taxon_phylogenetic_path_edges[desc1.taxon][desc2.taxon] = pedges
                    del (c1.desc_paths)
        self._mirror_lookups()