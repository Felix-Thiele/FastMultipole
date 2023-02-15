

from tools import *

import numpy as np
from QuadTree import QuadTree



class BarnesHut:

    def __init__(self, bodies):
        self.qtree = QuadTree(bodies, bodies_per_node=1)
        self.quad_dict = {}
        print(self.qtree.get_close_bodies(('00', '00')))

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

    def calc_forces(self):
        self. calc_node_centres()

        def bodies_in_node(node):
            if self.qtree.nodes[node][QuadTree.IS_LEAF]:
                return [list(_) for _ in self.qtree.nodes[node][QuadTree.BODIES]]
            else:
                bodies = []
                for child in self.qtree.children_labels(node):
                    bodies += bodies_in_node(child)
            return bodies

        forces = {}
        leaves = [node for node in self.qtree.nodes if self.qtree.nodes[node][QuadTree.IS_LEAF]]
        for leaf in leaves:
            leaf_details = self.qtree.nodes[leaf]
            for body in leaf_details[QuadTree.BODIES]:
                print('===================')
                print(body)
                print('====')
                force = 0
                # Potentials from bodies in leaf and neighboured nodes
                close_bodies = self.qtree.get_close_bodies(leaf)
                for close_body in close_bodies:
                    if not np.array_equal(close_body[:2], body[:2]):
                        print(close_body)
                        dist = complex(*body[:2]) - complex(*close_body[:2])
                        force += np.conjugate(close_body[2] * body[2] * 1 / dist)
                forces[tuple(body)] = force
                temp_leaf = leaf

                print('====')
                while len(temp_leaf[0])>1:
                    temp_leaf = temp_leaf[0][:-1], temp_leaf[1][:-1]
                    print(temp_leaf)
                    for neighbour in self.qtree.nodes[temp_leaf][QuadTree.ACT_NEIGH]:
                        close_bodies = bodies_in_node(neighbour)
                        for close_body in close_bodies:
                            if not np.array_equal(close_body[:2], body[:2]):
                                print(close_body)
                                dist = complex(*body[:2]) - complex(*close_body[:2])
                                force += np.conjugate(close_body[2] * body[2] * 1 / dist)
                        #dist = complex(*body[:2]) - complex(*self.quad_dict[neighbour][:2])
                        #force += np.conjugate(close_body[2] * body[2] * 1 / dist)
        return forces




