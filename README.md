UE computation, latency inference, toll pricing in traffic assignment
==========================


Setup
-----
Python dependencies (once only):

    sudo easy_install pip
    pip install cvxopt
    pip install networkx

Installing GDAL on a Mac: http://www.gis.usu.edu/~chrisg/python/2009/docs/gdal_mac.pdf

Mac OS X frameworks for GDAL: http://www.kyngchaos.com/software:frameworks

Running
-----
python draw_test.py

python read_shps_csv.py

python test_ue_solver.py

python shortest_paths_test.py

python inverse_opt_test.py

python isttt2014_experiments.py

python ecc2015_experiments.py

python toll_pricing_test.py

Troubleshooting
--------
...

Network of L.A.
--------

Coordinates for bounding box in L.A.: [-118.328299, 33.984601, -117.68132, 34.255881]

<img src="figures/map_graph.jpg" width=300px />

Add flow in equilibrium to model morning congestion

<img src="figures/map_congestion.jpg" width=300px />

References
--------

Patriksson's book: The Traffic Assignment Problem - Models and Methods

chapter 3 of urban networks: http://web.mit.edu/sheffi/www/selectedMedia/sheffi_urban_trans_networks.pdf

Anil’s paper: http://www.nt.ntnu.no/users/skoge/prost/proceedings/acc11/data/papers/1021.pdf

Toll pricing: http://www.columbia.edu/~ns2224/GT/hearnRamanaTollPricing.pdf

"A Multi-Convex approach to Latency Inference and Control in Traffic Equilibria from Sparse data" in preparation

"Approximate Bilevel Programming via Pareto Optimization for Imputation and Control of Optimization and Equilibrium models" in preparation
