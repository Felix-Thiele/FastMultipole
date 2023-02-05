# Some Fast Multipole explanations and implementations:

# [1] https://www.acom.rwth-aachen.de/_media/3teaching/0classes/ceswseminar/084_ces_seminar_siccha.pdf
# [2] https://www.math.uci.edu/~chenlong/226/FMMsimple.pdf
# [4] https://github.com/barbagroup/FMM_tutorial
# [3] https://github.com/lbluque/fmm
# [5] https://www.youtube.com/watch?v=qMLIyZi8Scenter_shift


from tools import *

import numpy as np
import itertools
import math
from scipy.special import binom

from QuadTree import QuadTree


# In the Fast Multipole Algorithm we want to approximate the forces between n bodies. The main idea is that we can
# approximate the forces between a node and some far away cluster of nodes. To define the clusters we use the QuadTree
# structure, a tree like clustur strukture. This defines a neighbour relation between cluster, and to approximate the
# force on a body we only need to calculate the forces from bodies in the same cluster, and then approximate the forces
# of neighbours of parent clusters. This algorithm is called barnes-hut algorithm and runs in O(n log(n)), but when
# implemented carefully, with the right approximation of forces between clusters can run in o(n) time and is then
# called the fast multipole algorithm. The Idea is to calculate the taylor approximations of forces outgoing and incoming
# into each cluster.
# The Outgoing coefficients can be calculated from the leaf upwords to the root, since the approximation
# of a parent can be calculated from the taylor approximations of its children clusters, by summing them up after shifting
# their origins.
# Outgoing coefficients can be recalulated as incoming coefficients of their neighbours in the same way, which again get
# summed up for each cluster.
# Finaly the part of the incoming approximations of a child that are due from its parents active neighbours can be
# added to the childs incoming approximation, by shifting the origin of the incoming taylor approximation to the center
# of the child. This is done recursively from the root to the child
# After all this we can calculate the force of each node by inserting its position into the inxoming taylor approx of the
# leaf cluster in wich it is contained.



class FastMultipole:

    def __init__(self, bodies):
        self.qtree = QuadTree(bodies, bodies_per_node=1)
        self.terms = 5
        self.outgoing_coeff_dict = {}
        self.incoming_coeff_dict = {}


    def _calc_outgoing_coef(self, node = ('', '')):
        # Recursive function to calculate the taylor approximations of forces outgoign from a quadrant around its center
        # For leaf nodes we dirrectly calculate them, and for parent nodes we sum up the taylor approximations
        #   of the children after shifting their centers to the center of the parent
        node_details = self.qtree.nodes[node]
        center = self.qtree.get_quad_center(node, node_details[QuadTree.LEVEL])
        if node_details[QuadTree.IS_LEAF]:
            self.outgoing_coeff_dict[node] = np.array([sum(b[2] for b in node_details[QuadTree.BODIES])] + [sum([-b[2]*complex(b[0] - center[0], b[1] - center[1])**k/k
                  for b in node_details[QuadTree.BODIES]]) for k in range(1, self.terms+1)])
        else:
            coefs = np.zeros((self.terms + 1), dtype=complex)
            for child in self.qtree.get_children(node):
                child_coeffs = self._calc_outgoing_coef(child)
                center_shift = complex(*self.qtree.get_quad_center(node, self.qtree.nodes[child][QuadTree.LEVEL])) - complex(*center)
                shifted_child_coeff = np.array([child_coeffs[0]] + [sum([child_coeffs[k] * center_shift ** (l - k) * binom(l - 1, k - 1) - (child_coeffs[0] * center_shift ** l) / l
                                  for k in range(1, l)]) for l in range(1, len(child_coeffs))])
                coefs += shifted_child_coeff
            self.outgoing_coeff_dict[node] = coefs
        return self.outgoing_coeff_dict[node]

    def _calc_incoming_coef_from_neigh_outgoing_coef(self, node, neighbour):
        # here we calculate the incoming coef due to the outgoing coefficients of an active neighbour
        nieghbour_coeffs = self.outgoing_coeff_dict[neighbour]
        center_shift = complex(*self.qtree.get_quad_center(node, self.qtree.nodes[node][QuadTree.LEVEL])) \
             - complex(*self.qtree.get_quad_center(neighbour, self.qtree.nodes[neighbour][QuadTree.LEVEL]))
        return np.array([sum([(nieghbour_coeffs[k] / center_shift ** k) * (-1) ** k for k in range(1, len(nieghbour_coeffs))]) +
                    nieghbour_coeffs[0] * np.log(-center_shift)] + [(1 / center_shift ** l) * sum([(nieghbour_coeffs[k] / center_shift ** k) * binom(l + k - 1, k - 1) * (-1) ** k
                                          for k in range(1, len(nieghbour_coeffs))]) - nieghbour_coeffs[0] / ((center_shift ** l) * l)
                     for l in range(1, len(nieghbour_coeffs))])

    def _calc_all_incoming_from_outgoing(self, node = ('', '')):
        # Recursive function recalculating all outgoing coefficients into incoming coefficients
        inc_coeff = np.zeros((self.terms + 1), dtype=complex)
        for index, act_neigh in enumerate(self.qtree.nodes[node][QuadTree.ACT_NEIGH]):
            inc_coeff += self._calc_incoming_coef_from_neigh_outgoing_coef(node, act_neigh)
        self.incoming_coeff_dict[node] = inc_coeff
        if not self.qtree.nodes[node][QuadTree.IS_LEAF]:
            for child in self.qtree.get_children(node):
                self._calc_all_incoming_from_outgoing(child)

    def _add_incoming_coef_to_children(self, node=('', '')):
        # Recursive function to calculate the incoming coeffients of a child due to its parents
        leaves = []
        node_details = self.qtree.nodes[node]
        if not node_details[QuadTree.IS_LEAF]:
            for child in self.qtree.get_children(node):
                center_shift = complex(*self.qtree.get_quad_center(node, self.qtree.nodes[node][QuadTree.LEVEL])) \
                               - complex(*self.qtree.get_quad_center(child, self.qtree.nodes[child][QuadTree.LEVEL]))
                node_inc_coeff = self.incoming_coeff_dict[node]
                shift_for_child_coeff = np.array([sum([node_inc_coeff[k]*binom(k,l)*(-center_shift)**(k-l)
                                          for k in range(l,len(node_inc_coeff))])
                                          for l in range(len(node_inc_coeff))])
                if child in self.incoming_coeff_dict:
                    self.incoming_coeff_dict[child] += shift_for_child_coeff
                else:
                    self.incoming_coeff_dict[child] = shift_for_child_coeff
                leaves += self._add_incoming_coef_to_children(child)
        else:
            return [node]
        return leaves

    def _calc_forces_for_bodies_in_leaves(self, leaves):
        forces = {}
        for leaf in leaves:
            leaf_details = self.qtree.nodes[leaf]
            z0, coeffs = complex(*self.qtree.get_quad_center(leaf, self.qtree.nodes[leaf][QuadTree.LEVEL])), self.incoming_coeff_dict[leaf]
            for body in leaf_details[QuadTree.BODIES]:
                force=0
                # Forces from innertayler expansions
                z = complex(*body[:2])
                force -= np.real(np.polyval(coeffs[::-1], z - z0))
                # Forces from bodies in leaf and neighboured nodes
                close_bodies = self.qtree.get_close_bodies(leaf)
                for close_body in close_bodies:
                    if not np.array_equal(close_body, body):
                        dist = np.sqrt((body[0] - close_body[0])**2 + (body[1] - close_body[1])**2)
                        force -= body[2]*np.log(dist)
                forces[tuple(body)]=force
        return forces

    def calc_forces(self):
        self._calc_outgoing_coef(node=('', ''))
        self._calc_all_incoming_from_outgoing(node=('', ''))
        leaves = self._add_incoming_coef_to_children(node=('', ''))
        forces = self._calc_forces_for_bodies_in_leaves(leaves)
        return forces


FM = FastMultipole(np.array([[1,2,3], [2,1,4], [1.2,2,3], [2,1.3,4], [1,2.6,3], [2.8,1,4]]))
forces = FM.calc_forces()