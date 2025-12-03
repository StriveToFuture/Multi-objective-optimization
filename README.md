An adaptive, reference-vector-based multi-objective optimizer derived from NSGA-II
1. What it is
DLF-RVEA (Dynamic Leader-Follower RVEA) is a lightweight yet powerful evolutionary algorithm for many-objective problems.
It keeps NSGA-II’s fast non-dominated sorting but replaces the crowding-distance stage with a reference-vector-guided dynamic leader–follower strategy, delivering a crisp, well-distributed Pareto front with low CPU cost.
2. Key idea in one sentence
Use reference vectors to on-the-fly assign leaders and followers, automatically balancing convergence toward the true PF and diversity across the entire front—even when objectives scale to 5, 10, 15+.
