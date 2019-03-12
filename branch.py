import os
import operator
from timeit import default_timer as timer
import networkx as nx
import common as c

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

    # initialize anti-graph
    ag = nx.create_empty_copy(g)

    # start timer for branching
    init = timer()
    TIMEOUT = 5 * 60

    # track size of best solution encountered
    best_sol = len(g)

    def bnb(g, ag, current=0, ub=-1, is_split=False):
        nonlocal best_sol

        # check if branching has timed out
        # if it has, error out with size of best solution encountered
        if timer() - init > TIMEOUT:
            raise ValueError(best_sol + len(slv))

        # check if graph is empty or a single vertex
        n = len(g)
        if n <= 1:
            return set()

        # check if graph is a single edge
        if n == 2 and g.size() == 1:
            return set([list(g)[0]])

        vc = set()

        # split instance into connected components if relevant
        if not nx.is_connected(g):
            # ccs = [(g.subgraph(cc).copy(), ag.subgraph(cc).copy())
            #        for cc in nx.connected_components(g)]

            ccs = [[v for v in g.subgraph(cc)] for cc in nx.connected_components(g)]

            # for (cc, acc) in ccs:
            #     vc = vc.union(bnb(cc, acc, current=current, ub=ub, is_split=True))

            for cc in ccs:
                excluded_vertices = list(set([v for v in g]) - set(cc))
                excluded_edges = []
                for v in excluded_vertices:
                    v_edges = [e for e in g.edges(v)]
                    excluded_edges += v_edges
                    for e in v_edges:
                        ag.add_edge(*e)
                    g.remove_node(v)

                vc = vc.union(bnb(g, ag, current=current, ub=ub, is_split=True))

                for v in excluded_vertices:
                    g.add_node(v)
                for e in excluded_edges:
                    g.add_edge(*e)
                    ag.remove_edge(*e)


            if not is_split and len(vc) < best_sol:
                best_sol = len(vc)

            return vc

        # get relevant instance info
        deg_list = sorted([(v, d) for v, d in g.degree()],
                          key=operator.itemgetter(1), reverse=True)
        v, max_deg = deg_list[0]
        _, min_deg = deg_list[-1]

        # apply degree-one reduction if relevant
        has_deg_one_redux = False
        if min_deg == 1:
            has_deg_one_redux = True
            tree_vc, cover = c.deg_one_redux(g, ag)
            vc = vc.union(tree_vc)
            current += len(tree_vc)

            if g.size() == 0:
                for e in cover:
                    g.add_edge(*e)
                    ag.remove_edge(*e)
                if not is_split and len(vc) < best_sol:
                    best_sol = len(vc)
                return vc

            deg_list = sorted([(v, d) for v, d in g.degree()],
                              key=operator.itemgetter(1), reverse=True)
            v, max_deg = deg_list[0]

        # simple bound check
        if ub > 0 and current >= ub:
            # add all vertices to the VC, effectively terminating the branch
            for v in g:
                vc.add(v)

            if has_deg_one_redux:
                for e in cover:
                    g.add_edge(*e)
                    ag.remove_edge(*e)

            return vc

        # bound check from Lemma 2.3 of Cygan et al.
        if max_deg <= ub and (len(g) > (ub**2 + ub) or g.size() > ub**2):
            # add all vertices to the VC, effectively terminating the branch
            for v in g:
                vc.add(v)

            if has_deg_one_redux:
                for e in cover:
                    g.add_edge(*e)
                    ag.remove_edge(*e)

            return vc

        # branch 1: v in VC
        in_cover = [e for e in g.edges(v)]
        in_vc = set([v])
        in_current = current + 1
        for e in in_cover:
            g.remove_edge(*e)
            ag.add_edge(*e)
        in_result = bnb(g, ag, current=in_current, ub=ub, is_split=is_split)
        in_current += len(in_result)
        in_vc = in_vc.union(in_result)
        if ub < 0 or in_current < ub:
            ub = in_current
        for e in in_cover:
            g.add_edge(*e)
            ag.remove_edge(*e)
        if not is_split and len(vc.union(in_vc)) < best_sol:
            best_sol = len(vc.union(in_vc))

        # branch 2: v is not in VC
        out_cover = []
        out_vc = set()
        out_current = current
        neighbors = [u for u in g[v]]
        for u in neighbors:
            out_subcover = [e for e in g.edges(u)]
            for e in out_subcover:
                out_cover.append(e)
                g.remove_edge(*e)
                ag.add_edge(*e)
            out_vc.add(u)
            out_current += 1
        out_result = bnb(g, ag, current=out_current, ub=ub, is_split=is_split)
        out_current += len(out_result)
        out_vc = out_vc.union(out_result)
        for e in out_cover:
            g.add_edge(*e)
            ag.remove_edge(*e)

        if has_deg_one_redux:
            for e in cover:
                g.add_edge(*e)
                ag.remove_edge(*e)

        if in_current < out_current:
            return vc.union(in_vc)
        else:
            vc = vc.union(out_vc)
            if not is_split and len(vc) < best_sol:
                best_sol = len(vc)
            return vc

    return slv.union(bnb(g, ag))


if __name__ == "__main__":
    with open("bnb_output.txt", 'w') as output:
        output.write("file_name vc_size run_time\n")
        with os.scandir(input_dir) as dir:
            for file in dir:
    # file = os.path.join(input_dir, "vc-exact_009.hgr")
                with open(file, encoding="latin-1") as f:
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
                    data = ' '.join([str(file.name), str(sol), str(run_time)])
                    output.write(data + '\n')
                    print(str(file.name) + " done")
                    print(sol, run_time)
