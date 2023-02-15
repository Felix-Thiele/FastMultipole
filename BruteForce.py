import numpy as np

class BruteForce:

    def __init__(self, bodies, force_fkt=lambda dist, q1, q2: np.conjugate(q1 * q2 * 1 / dist)):
        self.bodies = bodies
        self.force_fkt = force_fkt

    def calc_forces(self, force_fkt=None):
        if not force_fkt: force_fkt=self.force_fkt
        forces = {}
        for body in self.bodies:
            force = 0
            for sec_body in self.bodies:
                if not np.array_equal(sec_body[:2], body[:2]):
                    dist = complex(*body[:2]) - complex(*sec_body[:2])
                    force += force_fkt(dist, body[2], sec_body[2])
            forces[tuple(body)] = force
        return forces

    def calc_edge_forces(self, edges, force_fkt=None, edge_weight_influence=False):
        # edges are lists containing the two bodies and if edge_weight_influence also a weight factor
        if not force_fkt: force_fkt=self.force_fkt
        forces = {_: 0 for _ in self.bodies}
        for edge in edges:
            assert not np.array_equal(edge[0][:2], edge[1][:2])
            dist = complex(*edge[0][:2]) - complex(*edge[1][:2])
            weight_fac=edge[3] if edge_weight_influence else 1
            forces[tuple(edge[0])] += force_fkt(dist, edge[0][2], edge[1][2]) * weight_fac
            forces[tuple(edge[1])] -= force_fkt(dist, edge[0][2], edge[1][2]) * weight_fac
        return forces

    def gravity(self, gravity_factor=1, force_fkt=None):
        if not force_fkt: force_fkt=self.force_fkt
        center_gravity = np.mean(self.bodies[:2], axis=0)
        forces = {}
        for body in self.bodies:
            dist = complex(*body[:2]) - complex(*center_gravity)
            forces[tuple(body)] = force_fkt(dist, gravity_factor, body[2])