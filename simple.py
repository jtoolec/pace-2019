import os
import networkx as nx
import networkx.algorithms.approximation as naa
import common as c

script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")

# simple degree 1 reduction
def deg_one_redux(graph):
    vc = set()
    keep_going = True
    while keep_going:
        keep_going = False
        for v, d in graph.degree():
            if d == 1:
                u = [u for u in graph[v]][0]
                vc.add(u)
                leaves = []
                for w in graph[u]:
                    if graph.degree[w] == 1:
                        leaves.append(w)
                graph.remove_node(u)
                for w in leaves:
                    graph.remove_node(w)
                keep_going = True
                break
    return (vc, graph)

if __name__ == "__main__":
    with open("simple_output.txt", 'w') as output:
        output.write("file_name num_vertices num_edges " \
                   + "min_deg max_deg treewidth approx_vc\n")
        with os.scandir(input_dir) as dir:
            for file in dir:
                with open(file, encoding="latin-1") as f:
                    g = c.parse_graph(f)
                    vc, g = deg_one_redux(g)
                    if g.number_of_nodes() > 0:
                        v = g.number_of_nodes()
                        e = g.number_of_edges()
                        min_deg = min(d for v, d in g.degree())
                        max_deg = max(d for v, d in g.degree())
                        tw = naa.treewidth_min_degree(g)[0]
                        approx_vc = len(naa.min_weighted_vertex_cover(g))
                        data = ' '.join([str(file.name), str(v), str(e), \
                                         str(min_deg), str(max_deg), \
                                         str(tw), str(approx_vc)])
                        output.write(data + '\n')
                    else:
                        print("Vertex cover found for " + str(file.name))
                        print(len(vc))
