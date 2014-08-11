'''
Created on Aug 10, 2014

@author: jeromethai
'''


"""The goal is to compare different path solvers
Generate synthetic data:
1. compute the UE link flow using node-link flow formulation
2. compute the 3-shortest paths for each OD pair with weights the delays in UE (optimal)
3. compute the UE path flow using link-path formulation with the paths enumerated above
4. check that the resulting link flow vectors are the same (yes!)
5. generate a trajectory of waypoints for each used path by taking the closest waypoints to each link along the path
6. gathering routes with same waypoints sequence and sum up the flow for these routes

We suppose that we have:
1. Topology of the network, #lanes, #ff_delays, OD pairs, OD flows (OSM, google maps, census data, tomogravity)
2. Network is in steady-state equilibrium (UE) (common assumption in traffic assignment)
3. We know the latency functions (inverse equilibrium)
4. We know exactly the used routes (google maps, K-shortest paths on UE delays)
5. We have waypoint trajectories and flows and the set of routes associated to each
6. We have partial link flows from static sensors

How can we estimate the path flows? 4 experiments:
1. min ||P*x-l|| s.t. U*x=r, x>=0 with U OD-path incidence (most under-determined)
2. Compute link-flow UE, then compute feasible path flows (very under-determined)
3. Compute UE path flows directly (very under-determined)
4. min ||P*x-l|| s.t. U*x=r, x>=0 with U WP-path incidence (weakly to non under-determined)

Remarks:
1. UE doesn't enable to determine path flows
2. Turning ratios not accessible
3. Hence path flows is not really accessible (few route estimation in the literature)
4. That's why 'route flow estimation' is not used in the literature
"""


def synthetic_data(demand):
    pass


def experiment(demand):
    """One experiment"""


def main():
    pass


if __name__ == '__main__':
    pass