import ipdb
import random
import numpy as np
from scipy.sparse import csr_matrix

# Helper functions
# -------------------------------------
def to_np(X):
    return np.array(X).squeeze()

def to_sp(X):
    return csr_matrix((to_np(X.V),(to_np(X.I),to_np(X.J))), shape=X.size)

class LinkPath:
    def __init__(self, G, x_true, N=10):
        self.G = G
        self.x_true = x_true

        self.sample_linkpath(N=N)

    def sample_linkpath(self, N=10):
        self.lp = random.sample(self.G.links.keys(),N)
        self._get_lp_trajs()

    def _get_lp_trajs(self):
        rs = self.G.paths
        path_lps = [(r,[e.repr() for e in rs[r].links if e.repr() in self.lp]) \
                    for r in rs.keys()]
        lps = {}
        for value,key in path_lps:
            lps.setdefault(tuple(key), []).append(value)
        if () in lps:
            del lps[()]
        self.path_lps, self.lp_trajs = path_lps, lps

    def update_lp_flows(self):
        self.lp_flows = [sum([self.G.paths[i].flow for i in paths]) for \
                         paths in self.lp_trajs.values()]

    def simplex_lp(self):
        """Build simplex constraints from lp flows
        """
        from cvxopt import matrix, spmatrix
        n = len(self.lp_trajs)
        m = len(self.G.paths)
        if n == 0:
            return None, None
        I, J, r = [], [], matrix(0.0, (n,1))
        for i, path_ids in enumerate(self.lp_trajs.itervalues()):
            r[i] = self.lp_flows[i]
            for id in path_ids:
                I.append(i)
                J.append(self.G.indpaths[id])
        V = to_sp(spmatrix(1.0, I, J, (n, m)))
        r = to_np(r)
        return V, r

