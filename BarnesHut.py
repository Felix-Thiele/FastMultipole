

from tools import *

import numpy as np
from QuadTree import QuadTree



class BarnesHut:

    def __init__(self, bodies, force_fkt=lambda dist, q1, q2: np.conjugate(q1 * q2 * 1 / dist)):
        self.force_fkt = force_fkt
        self.qtree = QuadTree(bodies, bodies_per_node=1)
        self.quad_dict = {}

    def calc_node_centres(self, node = ('', '')):
        center, charge = 0, 0
        node_details = self.qtree.nodes[node]
        if node_details[QuadTree.IS_LEAF]:
            for b in node_details[QuadTree.BODIES]:
                center += complex(b[0], b[1])
                charge += b[2]
            self.quad_dict[node] = (center, charge)
        else:
            for child in self.qtree.children_labels(node):
                ce, ch = self.calc_node_centres(child)
                center += ce
                charge += ch
        self.quad_dict[node]=(center, charge)
        return center, charge

    def calc_forces(self, force_fkt=None):

        if not force_fkt: force_fkt = self.force_fkt

        self. calc_node_centres()

        forces = {}
        leaves = [node for node in self.qtree.nodes if self.qtree.nodes[node][QuadTree.IS_LEAF]]
        for leaf in leaves:
            leaf_details = self.qtree.nodes[leaf]
            for body in leaf_details[QuadTree.BODIES]:
                force = 0
                # Potentials from bodies in leaf and neighboured nodes
                close_bodies = self.qtree.get_close_bodies(leaf)
                for close_body in close_bodies:
                    if not np.array_equal(close_body[:2], body[:2]):
                        dist = complex(*body[:2]) - complex(*close_body[:2])
                        force += np.conjugate(close_body[2] * body[2] * 1 / dist)
                temp_leaf = leaf
                while len(temp_leaf[0])>1:
                    for neighbour in self.qtree.nodes[temp_leaf][QuadTree.ACT_NEIGH]:
                        # This is the same as brute force to test QTREE
                        """close_bodies = self.qtree.bodies_in_node(neighbour)
                        for close_body in close_bodies:
                            if not np.array_equal(close_body[:2], body[:2]):
                                dist = complex(*body[:2]) - complex(*close_body[:2])
                                force += force_fkt(dist, body[2], close_body[2])"""
                        close_body, charge = self.quad_dict[neighbour]
                        dist = complex(*body[:2]) - close_body
                        force += force_fkt(dist, body[2], charge)

                    temp_leaf = temp_leaf[0][:-1], temp_leaf[1][:-1]

                forces[tuple(body)] = force
        return forces




