import os
import operator
import networkx as nx
import common as c

script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")

# wrapper function to keep track of best VC found during bnb
def solve(g):
    ub = -1
    slv = set()

    # include vertices with self-loops in any VC
    for e in nx.selfloop_edges(g):
        v, _ = e
        slv.add(v)

    # self-loop vertices can be removed from the graph
    for v in slv:
        g.remove_node(v)

    def bnb(g, current=0):
        nonlocal ub

        # check if graph is empty or a single vertex
        n = len(g)
        if n <= 1:
            if ub < 0 or current < ub:
                ub = current
            return set()

        # check if graph is a single edge
        if n == 2 and g.size() == 1:
            current += 1
            if ub < 0 or current < ub:
                ub = current
            return set([list(g)[0]])

        vc = set()

        # split instance into connected components if relevant
        if not nx.is_connected(g):
            ccs = [g.subgraph(cc) for cc in nx.connected_components(g)]

            for cc in ccs:
                vc = vc.union(bnb(cc, current=current))

            return vc

        # get relevant instance info
        deg_list = sorted([(v, d) for v, d in g.degree()],
                          key=operator.itemgetter(1), reverse=True)
        v, max_deg = deg_list[0]
        _, min_deg = deg_list[-1]

        # apply degree-one reduction if relevant
        if min_deg == 1:
            g, tree_vc = c.deg_one_redux(g)
            vc = vc.union(tree_vc)
            current += len(tree_vc)

            if len(g) == 0:
                if ub < 0 or current < ub:
                    ub = current
                return vc

            deg_list = sorted([(v, d) for v, d in g.degree()],
                              key=operator.itemgetter(1), reverse=True)
            v, max_deg = deg_list[0]

        # simple bound check
        if ub > 0 and current >= ub:
            # add all vertices to the VC, effectively terminating the branch
            for v in g:
                vc.add(v)

            return vc

        # branch 1: v in VC
        in_vc = vc.copy()
        in_g = g.copy()
        in_vc.add(v)
        in_g.remove_node(v)
        in_current = current + 1

        # branch 2: v is not in VC
        out_vc = vc.copy()
        out_g = g.copy()
        out_current = current
        for u in g[v]:
            out_vc.add(u)
            out_g.remove_node(u)
            out_current += 1
        out_g.remove_node(v)

        # recursion
        in_result = in_vc.union(bnb(in_g, current=in_current))
        out_result = out_vc.union(bnb(out_g, current=out_current))
        if len(in_result) < len(out_result):
            return in_result
        else:
            return out_result

    return slv.union(bnb(g))


if __name__ == "__main__":
    # with open("bnb_output.txt", 'w') as output:
    #     output.write("file_name vc_size\n")
    #     with os.scandir(input_dir) as dir:
    #         for file in dir:
    file = os.path.join(input_dir, "vc-exact_199.hgr")
    with open(file, encoding="latin-1") as f:
        g = c.parse_graph(f)
        vc = solve(g)
        # data = ' '.join([str(file.name), len(vc)])
        # output.write(data + '\n')
        # print(str(file.name) + " done")
        print(len(vc))
