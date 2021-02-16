import networkx as nx
from webweb import Web

G = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')
HPO = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')

# remove connection between CFTR and it's known HPO terms
terms_to_remove = ['HP:0011949', 'HP:0001733', 'HP:0002202', 'HP:0030830', 'HP:0000006', 'HP:0011962', 'HP:0002105',
                   'HP:0000007', 'HP:0004313', 'HP:0002027', 'HP:0002721', 'HP:0005952', 'HP:0006528', 'HP:0006538',
                   'HP:0002110', 'HP:0030877', 'HP:0002613', 'HP:0011961', 'HP:0001738', 'HP:0030247', 'HP:0002783',
                   'HP:0000118', 'HP:0000837', 'HP:0005376', 'HP:0001508', 'HP:0001658', 'HP:0003251', 'HP:0000819',
                   'HP:0012873', 'HP:0100749', 'HP:0002099', 'HP:0005213', 'HP:0002024', 'HP:0004401', 'HP:0000798',
                   'HP:0001425', 'HP:0002205', 'HP:0100027', 'HP:0012379', 'HP:0001945', 'HP:0004469', 'HP:0011947',
                   'HP:0002240', 'HP:0002150', 'HP:0011227', 'HP:0001974', 'HP:0004326', 'HP:0002570', 'HP:0001648',
                   'HP:0002035']

for hpo in terms_to_remove:
    try:
        G.remove_edge(hpo, 'CFTR')
    except nx.exception.NetworkXError:
        pass

paths = []
# length of paths getting from those terms to CFTR now
for hpo in terms_to_remove:
    try:
        temp_paths = nx.all_shortest_paths(G, hpo, 'CFTR')
    except:
        print('No path from ' + hpo + ' to CFTR')
        continue
    if len(list(temp_paths)[0]) > 3:
        print(hpo)
    # print(temp_paths)
    for path in temp_paths:
        # print(len(path))
        # show the visualization
        paths.append(path)

nodes = []
for path in paths:
    for p in path:
        nodes.append(p)

g = nx.subgraph(G, nodes)

display = {
    'nodes': {str(i): {'name': node, 'deg': G.degree(node), 'related': node in terms_to_remove} for i, node in enumerate(g.nodes())}
}

web = Web(nx.to_numpy_array(g), display=display)
web.display.colorBy = 'related'
web.display.sizeBy = 'deg'
web.display.showNodeNames = True
web.display.linkLength = 100
web.display.charge = 300
web.display.radius = 10
web.show()

# get depth down HPO
depth = []
for hpo in terms_to_remove:
    depth.append(len(nx.shortest_path(HPO,hpo,'HP:0000001')))
    nx.degree(g,hpo)