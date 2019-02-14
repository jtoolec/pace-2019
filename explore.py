import os
import networkx as nx
import networkx.algorithms.approximation as naa

# get basic graph info listed in file header
script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")
with open("basic.txt", 'w') as output:
    output.write("file_name num_vertices num_edges avg_deg\n")
    with os.scandir(input_dir) as dir:
        for file in dir:
            with open(file, encoding="latin-1") as f:
                for line in f:
                    if line[0] == 'p':
                        s = line.split()
                        v = s[2]
                        e = s[3]
                        avg_deg = int(s[3]) / int(s[2])
                        data = ' '.join([str(file.name), str(v), str(e), str(avg_deg)])
                        output.write(data + '\n')
                        break

# get more detailed information on problem instances
with open("detail.txt", 'w') as output:
    output.write("file_name min_deg max_deg treewidth approx_vc\n")
    with os.scandir(input_dir) as dir:
        for file in dir:
            g = nx.Graph()
            with open(file, encoding="latin-1") as f:
                for line in f:
                    if line[0] == 'p':
                        s = line.split()
                        v = int(s[2])
                        g.add_nodes_from(range(1, v + 1))
                    elif line[0] != 'c':
                        e = line.split()
                        g.add_edge(int(e[0]), int(e[1]))
                min_deg = min(d for n, d in g.degree())
                max_deg = max(d for n, d in g.degree())
                tw = naa.treewidth_min_degree(g)[0]
                vc = len(naa.min_weighted_vertex_cover(g))
                data = ' '.join([str(file.name), str(min_deg), str(max_deg), str(tw), str(vc)])
                output.write(data + '\n')
