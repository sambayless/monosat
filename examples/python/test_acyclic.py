import networkx as nx
import networkx.algorithms.tree
import random
from monosat import *
from networkx.relabel import convert_node_labels_to_integers

print("begin encode");

seed = random.randint(1,100000)

random.seed(seed)
print("RandomSeed=" + str(seed))

size =6

#Monosat().newSolver("-decide-graph-bv -no-decide-theories -no-decide-graph-rnd   -lazy-maxflow-decisions -conflict-min-cut -conflict-min-cut-maxflow -reach-underapprox-cnf ")




graphs = [nx.diamond_graph(), nx.complete_graph(2), nx.complete_graph(3), nx.complete_graph(4), nx.cycle_graph(5),
          nx.grid_2d_graph(4,4),nx.grid_2d_graph(4,4,True),
          nx.dense_gnm_random_graph(4,8,123),
          nx.dense_gnm_random_graph(4,8,124),
          nx.dense_gnm_random_graph(4,16,125)
          ]
            #Directed graphs:
graphs+=[nx.gnp_random_graph(4,8,123,True),nx.gnp_random_graph(4,8,124,True),nx.gnp_random_graph(4,8,125,True),nx.gnp_random_graph(4,8,126,True),
         nx.gnp_random_graph(4,16,127,True),
          nx.gnp_random_graph(4,32,128,True),
          nx.gnp_random_graph(9,16,127,True),
         
         ]
for i,nxg in enumerate(graphs):

    nxg = convert_node_labels_to_integers(nxg)
    g = Graph()
    print("Graph %d"%(i))
    nodes = []
    for n in nxg:
        nodes.append(g.addNode())
    nedges=0
    for u,v in nxg.edges():

        Assert(g.addEdge(nodes[u],nodes[v]))

        print("%d -> %d"%(u,v))

        #if(nxg.is_directed()):
        #    Assert(g.addEdge(nodes[v],nodes[u]))
        
    if(nxg.is_directed()):
        if(nx.is_directed_acyclic_graph(nxg)):
            print("Acyclic")
            Assert(g.acyclic(True))
        else:
            print("Has cycle")
            Assert(Not(g.acyclic(True)))

    else:
        if(networkx.algorithms.tree.recognition.is_forest(nxg)):
            print("Undirected Acyclic")
           
            Assert(g.acyclic(False))
        else:
            print("Undirected Has cycle")
            
            Assert(Not(g.acyclic(False)))
        

for i,nxg in enumerate(graphs):

    nxg = convert_node_labels_to_integers(nxg)
    g = Graph()
    print("Graph %d"%(i))
    nodes = []
    for n in nxg:
        nodes.append(g.addNode())
    nedges=0
    for u,v in nxg.edges():

        g.addEdge(nodes[u],nodes[v])

        print("%d -> %d"%(u,v))

        #if(nxg.is_directed()):
        #    Assert(g.addEdge(nodes[v],nodes[u]))
        
    if(nxg.is_directed()):
        if(nx.is_directed_acyclic_graph(nxg)):
            print("Acyclic")
            Assert(g.acyclic(True))
        else:
            print("Has cycle")
            Assert(Not(g.acyclic(True)))

    else:
        if(networkx.algorithms.tree.recognition.is_forest(nxg)):
            print("Undirected Acyclic")
           
            Assert(g.acyclic(False))
        else:
            print("Undirected Has cycle")
            
            Assert(Not(g.acyclic(False)))
        


result = Solve()
print("Result is " + str(result))
assert(result==True)

