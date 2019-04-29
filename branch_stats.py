import os
import sys
import operator
from timeit import default_timer as timer
import networkx as nx
import networkx.algorithms.approximation as naa
import common as c

sys.setrecursionlimit(2000)

script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")

# wrapper function to keep track of best VC found during bnb
def solve(g):
    slv = set()

    # include vertices with self-loops in any VC
    for e in nx.selfloop_edges(g):
        v, _ = e
        slv.add(v)

    # self-loop vertices can be removed from the graph
    for v in slv:
        g.remove_node(v)

    # run degree-one reduction once if relevant
    # this call of the reduction removes vertices to handle large trees quickly
    for _, d in g.degree():
        if d == 1:
            tree_vc = c.deg_one_redux(g, in_place=False)
            slv = slv.union(tree_vc)

            if g.size() == 0:
                return slv

            break

    # get 2-approx estimate for VC
    approx_vc = naa.min_weighted_vertex_cover(g)
    approx_ub = len(approx_vc)

    # start timer for branching
    init = timer()
    TIMEOUT = 2 * 60

    # track size of best solution encountered
    best_sol = len(g)

    def bnb(g, current=0, ub=-1, layer=0, is_split=False):
        nonlocal best_sol
        global base_time, deg_one_time, deg_list_time, simple_bound_time, cc_time, cygan_time, mirror_time, satellite_time, mv_branch_time, nsv_branch_time, cleanup_time
        global deg_one_count, simple_bound_count, cygan_count, mirror_count, satellite_count

        # check if branching has timed out
        # if it has, error out with size of best solution encountered
        if timer() - init > TIMEOUT:
            raise ValueError(best_sol + len(slv))

        base_start = timer()

        # check if graph is empty or a single vertex
        n = len(g)
        if n <= 1:
            base_time += timer() - base_start
            return set()

        # check if graph is a single edge
        if n == 2 and g.size() == 1:
            base_time += timer() - base_start
            return set([list(g)[0]])

        base_time += timer() - base_start

        vc = set()

        deg_one_start = timer()

        # apply degree-one reduction if relevant
        has_deg_one_redux = False
        for _, d in g.degree():
            if d == 1:
                deg_one_count += 1

                has_deg_one_redux = True
                tree_vc, cover = c.deg_one_redux(g)
                vc = vc.union(tree_vc)
                current += len(tree_vc)

                if g.size() == 0:
                    for e in cover:
                        g.add_edge(*e)
                    if not is_split and len(vc) < best_sol:
                        best_sol = len(vc)
                    deg_one_time += timer() - deg_one_start
                    return vc

                break

        deg_one_time += timer() - deg_one_start

        deg_list_start = timer()

        # get relevant instance info
        deg_list = sorted([(v, d) for v, d in g.degree()],
                          key=operator.itemgetter(1), reverse=True)
        v, max_deg = deg_list[0]

        deg_list_time += timer() - deg_list_start

        simple_bound_start = timer()

        # simple bound check
        if ub > 0 and current >= ub:
            simple_bound_count += 1

            # add all vertices to the VC, effectively terminating the branch
            for v in g:
                vc.add(v)

            if has_deg_one_redux:
                for e in cover:
                    g.add_edge(*e)

            simple_bound_time += timer() - simple_bound_start
            return vc

        simple_bound_time += timer() - simple_bound_start

        cc_start = timer()

        # split instance into connected components if relevant
        if not nx.is_connected(g):
            # ccs = [(g.subgraph(cc).copy(), ag.subgraph(cc).copy())
            #        for cc in nx.connected_components(g)]

            ccs = [[v for v in g.subgraph(cc)] for cc in nx.connected_components(g)]

            # for (cc, acc) in ccs:
            #     vc = vc.union(bnb(cc, acc, current=current, ub=ub, is_split=True))

            for cc in ccs:
                if len(cc) > 1:
                    excluded_vertices = list(set([v for v in g]) - set(cc))
                    excluded_edges = []
                    for v in excluded_vertices:
                        v_edges = [e for e in g.edges(v)]
                        excluded_edges += v_edges
                        g.remove_node(v)

                    cc_time += timer() - cc_start

                    cc_vc = bnb(g, current=current, ub=-1, layer=layer, is_split=True)

                    cc_start = timer()

                    vc = vc.union(cc_vc)
                    current += len(cc_vc)

                    for v in excluded_vertices:
                        g.add_node(v)
                    for e in excluded_edges:
                        g.add_edge(*e)

            if has_deg_one_redux:
                for e in cover:
                    g.add_edge(*e)

            if not is_split and len(vc) < best_sol:
                best_sol = len(vc)

            cc_time += timer() - cc_start
            return vc

        cc_time += timer() - cc_start

        cygan_start = timer()

        # bound check from Lemma 2.3 of Cygan et al.
        if max_deg <= ub and (len(g) > (ub**2 + ub) or g.size() > ub**2):
            cygan_count += 1

            # add all vertices to the VC, effectively terminating the branch
            for v in g:
                vc.add(v)

            if has_deg_one_redux:
                for e in cover:
                    g.add_edge(*e)

            cygan_time += timer() - cygan_start
            return vc

        cygan_time += timer() - cygan_start

        mirror_start = timer()

        # find mirrors of v
        neighbors = [u for u in g[v]]
        second_neighbors = c.neighborhood(g, v, 2)
        mirrors = [v] # always include v itself for branching
        for u in second_neighbors:
            neighbor_diff = list(set(neighbors) - set([w for w in g[u]]))
            s = len(neighbor_diff)
            if s == 0:
                mirror_count += 1

                mirrors.append(u)
            else:
                max_edges = (s * (s + 1)) / 2
                if g.subgraph(neighbor_diff).size() == max_edges:
                    mirror_count += 1

                    mirrors.append(u)

        mirror_time += timer() - mirror_start

        satellite_start = timer()

        # if no mirrors, find satellites of v
        satellites = set([v])
        if len(mirrors) == 1:
            for w in neighbors:
                neighbor_diff = list(set([u for u in g[w]]) - set(neighbors + [v]))
                if len(neighbor_diff) == 1:
                    satellite_count += 1

                    satellites.add(neighbor_diff[0])

        satellite_time += timer() - satellite_start

        mv_branch_start = timer()

        # branch 1: M[v] in VC
        in_cover = []
        in_vc = set(mirrors)
        in_current = current + len(in_vc)
        for u in mirrors:
            in_subcover = [e for e in g.edges(u)]
            for e in in_subcover:
                in_cover.append(e)
                g.remove_edge(*e)

        mv_branch_time += timer() - mv_branch_start

        in_result = bnb(g, current=in_current, ub=ub, layer=layer+1, is_split=is_split)

        mv_branch_start = timer()

        in_current += len(in_result)
        in_vc = in_vc.union(in_result)
        if ub < 0 or in_current < ub:
            ub = in_current
        for e in in_cover:
            g.add_edge(*e)
        if not is_split and len(vc.union(in_vc)) < best_sol:
            best_sol = len(vc.union(in_vc))

        mv_branch_time += timer() - mv_branch_start

        nsv_branch_start = timer()

        # branch 2: N(S[v]) in VC
        out_cover = []
        out_vc = set()
        out_current = current
        satellite_neighbors = set()
        for u in satellites:
            for w in g[u]:
                satellite_neighbors.add(w)
        for w in satellite_neighbors:
            out_subcover = [e for e in g.edges(w)]
            for e in out_subcover:
                out_cover.append(e)
                g.remove_edge(*e)
            out_vc.add(w)
            out_current += 1

        nsv_branch_time += timer() - nsv_branch_start

        out_result = bnb(g, current=out_current, ub=ub, layer=layer+1, is_split=is_split)

        nsv_branch_start = timer()

        out_current += len(out_result)
        out_vc = out_vc.union(out_result)
        for e in out_cover:
            g.add_edge(*e)
        # if not is_split and len(vc.union(out_vc)) < best_sol:
        #     best_sol = len(vc.union(out_vc))

        nsv_branch_time += timer() - nsv_branch_start

        cleanup_start = timer()

        if has_deg_one_redux:
            for e in cover:
                g.add_edge(*e)

        if in_current < out_current:
            cleanup_time += timer() - cleanup_start
            return vc.union(in_vc)
        else:
            vc = vc.union(out_vc)
            if not is_split and len(vc) < best_sol:
                best_sol = len(vc)
            cleanup_time += timer() - cleanup_start
            return vc

    return slv.union(bnb(g, ub=approx_ub))


if __name__ == "__main__":
    with open("bnb_stats_output.txt", 'w') as output:
        output.write("file_name vc_size run_time base_time deg_one_time deg_list_time simple_bound_time " +
                     "cc_time cygan_time mirror_time satellite_time mv_branch_time nsv_branch_time " +
                     "cleanup_time deg_one_count simple_bound_count cygan_count mirror_count " +
                     "satellite_count\n")
        with os.scandir(input_dir) as dir:
            for file in dir:
    # file = os.path.join(input_dir, "vc-exact_197.hgr")
                with open(file, encoding="latin-1") as f:
                    base_time = 0
                    deg_one_time, deg_one_count = 0, 0
                    deg_list_time = 0
                    simple_bound_time, simple_bound_count = 0, 0
                    cc_time = 0
                    cygan_time, cygan_count = 0, 0
                    mirror_time, mirror_count = 0, 0
                    satellite_time, satellite_count = 0, 0
                    mv_branch_time = 0
                    nsv_branch_time = 0
                    cleanup_time = 0

                    start = timer()
                    g = c.parse_graph(f)
                    try:
                        vc = solve(g)
                        sol_found = True
                    except ValueError as err:
                        sol_found = False
                        vc = err.args[0]
                    except RuntimeError:
                        sol_found = False
                        vc = -1
                    end = timer()

                    if sol_found:
                        sol = len(vc)
                    else:
                        sol = vc
                    run_time = end - start
                    stats = [base_time, deg_one_time, deg_list_time, simple_bound_time,
                             cc_time, cygan_time, mirror_time, satellite_time,
                             mv_branch_time, nsv_branch_time, cleanup_time,
                             deg_one_count, simple_bound_count, cygan_count,
                             mirror_count, satellite_count]
                    stats = [str(round(stat, 5)) for stat in stats]
                    data = ' '.join([str(file.name), str(sol), str(run_time)] + stats)
                    output.write(data + '\n')
                    print(str(file.name) + " done")
                    print(sol, run_time, stats)
