import os
import networkx as nx

script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")

def parse_graph(file):
    g = nx.Graph()
    for line in file:
        if line[0] == 'p':
            s = line.split()
            v = int(s[2])
            g.add_nodes_from(range(1, v + 1))
        elif line[0] != 'c':
            e = line.split()
            g.add_edge(int(e[0]), int(e[1]))
    return g

def has_self_loop(file):
    for line in file:
        if line[0] != 'p' and line[0] != 'c':
            e = line.split()
            if int(e[0]) == int(e[1]):
                return True
    return False

# return nth neighborhood of a vertex
# taken from https://stackoverflow.com/questions/22742754/finding-the-n-degree-neighborhood-of-a-node
def neighborhood(g, v, n):
    path_lengths = nx.shortest_path_length(g, source=v)
    return [u for u, l in path_lengths.items() if l == n]

# simple degree 1 reduction
# also solves trees
def deg_one_redux(g, ag):
    # g = g.copy()
    vc = set()
    keep_going = True
    cover = []
    while keep_going:
        keep_going = False
        for v, d in g.degree():
            if d == 1:
                u = [u for u in g[v]][0]
                vc.add(u)
                # leaves = []
                # for w in g[u]:
                #     if g.degree[w] == 1:
                #         leaves.append(w)
                # g.remove_node(u)
                # for w in leaves:
                #     g.remove_node(w)
                subcover = [e for e in g.edges(u)]
                for e in subcover:
                    cover.append(e)
                    g.remove_edge(*e)
                    ag.add_edge(*e)
                keep_going = True
                break
    return (vc, cover)

if __name__ == "__main__":
    with os.scandir(input_dir) as dir:
        for file in dir:
            with open(file, encoding="latin-1") as f:
                if has_self_loop(f):
                    print(file.name)
