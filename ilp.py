import os
from timeit import default_timer as timer
import pulp

script_dir = os.path.dirname(__file__)
input_dir = os.path.join(script_dir, "vc-exact-public")

with open("ilp_output.txt", 'w') as output:
    output.write("file_name vc_size\n")
    with os.scandir(input_dir) as dir:
        for file in dir:
            with open(file, encoding="latin-1") as f:
                start = timer()
                prob = pulp.LpProblem("VC", pulp.LpMinimize)
                for line in f:
                    if line[0] == 'p':
                        s = line.split()
                        v = int(s[2])
                        vs = pulp.LpVariable.dicts("vertices", [n for n in range(1, v + 1)],
                                                   lowBound=0, upBound=1, cat="Integer")
                        # add objective function
                        prob += pulp.lpSum([vs[i] for i in range(1, v + 1)])
                    elif line[0] != 'c':
                        e = line.split()
                        # add edge constraint
                        prob += pulp.lpSum(vs[int(e[0])] + vs[int(e[1])]) >= 1
                #pulp.LpSolverDefault.msg = 1
                try: # attempt to solve ILP
                    prob.solve()
                    end = timer()
                    run_time = end - start
                    if prob.status != 1: # check if prob is optimal
                        data = ' '.join([str(file.name), '0'])
                    else:
                        data = ' '.join([str(file.name),
                                         str(int(pulp.value(prob.objective))),
                                         str(run_time)])
                    output.write(data + '\n')
                except: # continue even if solver errors out
                    data = ' '.join([str(file.name), '-1'])
                    output.write(data + '\n')
                print(str(file.name) + " done")
