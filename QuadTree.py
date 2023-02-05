from tools import *

import itertools
import numpy as np


#todo maybe save min and width of leaf individualy
#todo a bit sloppy the lvl is the length of the node, redundant
class QuadTree:
    LEVEL, IS_LEAF, ACT_NEIGH, BODIES = 0, 1, 2, 3

    # get the bodies as a n*3 np array with the x and y positions and the charge
    def __init__(self, bodies, bodies_per_node=1):
        assert len(np.unique(bodies, axis=0)) == len(bodies)
        # self.nodes is a dict containing all tree nodes, mapped to 
        #   - the level of the node
        #   - bool if it is a leaf
        #   - a list of interaction neighboutrs
        #   - the list of bodies it contains if it is a leaf
        self.nodes = {('', ''): [0, True, [], bodies]}
        self.x_min, self.y_min, self.width, self.height = np.min(bodies[:, 0]), np.min(bodies[:, 1]), \
            np.max(bodies[:, 0]) - np.min(bodies[:, 0]), np.max(bodies[:, 1]) - np.min(bodies[:, 1])

        def recursive_build(node, bodies_per_node=bodies_per_node):
            if len(self.nodes[node][QuadTree.BODIES]) > bodies_per_node:
                children = self._split_node(node)
                for child in children:
                    recursive_build(child)

        recursive_build(('', ''))

    def _get_interaction_nodes(self, node, lvl):
        # Get the interaction neighbors of a node
        interaction_nodes = []
        for x_shift, y_shift in itertools.product([-1, 0, 1], [-1, 0, 1]):
            x_shifted_parent, y_shifted_parent = bin_int(node[0][:-1]) - x_shift, bin_int(node[1][:-1]) - y_shift
            for x_child, y_child in itertools.product([0, 1], [0, 1]):
                if abs(2 * x_shifted_parent + x_child - bin_int(node[0])) > 1 or abs(
                        2 * y_shifted_parent + y_child - bin_int(node[1])) > 1:
                    x = int_bin(x_shifted_parent, lvl - 1) + str(x_child)
                    y = int_bin(y_shifted_parent, lvl - 1) + str(y_child)
                    interaction_nodes.append((x, y))
        return [n for n in interaction_nodes if n in self.nodes]

    def get_close_bodies(self, node):
        # Gets the direct neighbours of a node
        close_bodies = []
        for x_shift, y_shift in itertools.product([-1, 0, 1], [-1, 0, 1]):
            x_shifted, y_shifted = int_bin(bin_int(node[0])-x_shift, len(node[0])), int_bin(bin_int(node[1]) - y_shift, len(node[0]))
            while (x_shifted, y_shifted) not in self.nodes:
                x_shifted, y_shifted = x_shifted[:-1], y_shifted[:-1]
            close_bodies += [list(_) for _ in self.nodes[(x_shifted, y_shifted)][QuadTree.BODIES]]
        return np.unique(np.array(close_bodies), axis=0)

    def get_quad_center(self, node, lvl):
        return self.x_min + self.width *1 / 2 ** (lvl) * (1 / 2 + bin_int(node[0])), self.y_min + self.height * 1 / 2 ** (lvl) * (1 / 2 + bin_int(node[1]))

    def get_children(self, node):
        # This will return the codes of children independent of if they are already defined
        return [(node[0] + str(x), node[1] + str(y)) for x, y in itertools.product([0, 1], [0, 1])]

    def _split_node(self, node):
        new_nodes = []
        lvl = self.nodes[node][QuadTree.LEVEL]
        orig_bodies = self.nodes[node][QuadTree.BODIES]
        for new_node in self.get_children(node):
            bodies = orig_bodies.copy()
            new_nodes.append(new_node)
            # bodies of children
            x_split, y_split = self.get_quad_center(node, lvl)
            x_mask, y_mask = bodies[:, 0] < x_split, bodies[:,1] < y_split
            if new_node[0][-1]=='1': x_mask = np.logical_not(x_mask)
            if new_node[1][-1]=='1': y_mask = np.logical_not(y_mask)
            # interaction nodes of new child
            interaction_nodes = self._get_interaction_nodes(new_node, lvl + 1)
            self.nodes[new_node] = [lvl + 1,
                                    True,
                                    interaction_nodes,
                                    bodies[np.where(x_mask & y_mask)]]
            # add new child as interaction node to its interaction neighbours
            for interact in interaction_nodes:
                self.nodes[interact][QuadTree.ACT_NEIGH].append(new_node)
        self.nodes[node][QuadTree.IS_LEAF] = False
        self.nodes[node][QuadTree.BODIES] = []
        return new_nodes
