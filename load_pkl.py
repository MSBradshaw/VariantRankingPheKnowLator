import networkx as nx


def load_network(el, node_mapping=None):
    if node_mapping is not None:
        url_to_common_name = {line.strip().split()[0]: line.strip().split()[1] for line in open(node_mapping)}
    else:
        url_to_common_name = None
    G = nx.Graph()
    types = set()
    for line in open(el, 'r'):
        row = line.strip().split('\t')
        if len(row) == 3:
            # when using 'Edgelists/pkl_triples_with_symbols.tsv'
            G.add_edge(row[0], row[2])
            G.edges[row[0], row[2]]['type'] = row[1]
            G.edges[row[0], row[2]]['weight'] = .001
            types.add(row[0])
        elif len(row) == 4:
            subject = row[1]
            the_object = row[3]
            if subject in url_to_common_name:
                subject = url_to_common_name[subject]
            if the_object in url_to_common_name:
                the_object = url_to_common_name[the_object]
            # when using 'Edgelists/PKT_Master_Edge_List_NoOntologyData.txt'
            G.add_edge(subject, the_object)
            G.edges[subject, the_object]['type'] = row[2]
            G.edges[subject, the_object]['type_name'] = row[0]
            G.edges[subject, the_object]['weight'] = .001
            types.add(row[0])
    return G

# Load the ontology free version of PKL (only has expert curated edges)
G_ontology_free = load_network('Edgelists/PKT_Master_Edge_List_NoOntologyData.txt', 'Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt')

# load the original version of PKL, but with Gene common names instead of protein identifiers
G_original = load_network('Edgelists/pkl_triples_with_symbols.tsv')

