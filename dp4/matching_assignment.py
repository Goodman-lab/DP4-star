from collections import Counter
import networkx as nx

def external_scale_proton_shifts(calculated_shifts):

    scaled = [0.9770793502768845 * calc - 0.019505417520415236 for calc in calculated_shifts]

    return scaled

def list_difference(list1, list2):
    c1 = Counter(list1)
    c2 = Counter(list2)

    diff = c1-c2
    return list(diff.elements())

def find_assignment(exp, calc, threshold=0.6):

    scaled = external_scale_proton_shifts(calc)

    # Create a bipartite graph
    G = nx.Graph()

    # Add nodes for calc and exp with a bipartite attribute
    calc_nodes = [('calc', i) for i in range(len(scaled))]
    exp_nodes = [('exp', i) for i in range(len(exp))]
    G.add_nodes_from(calc_nodes, bipartite=0)
    G.add_nodes_from(exp_nodes, bipartite=1)

    # Add edges for all pairs within the threshold
    for i, c in enumerate(scaled):
        for j, e in enumerate(exp):
            if abs(c - e) <= threshold:
                G.add_edge(('calc', i), ('exp', j), weight=-abs(c - e))

    # Find the maximum matching
    matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

    matched_pair_indices = set()

    for pair in matching:
        # pair could be (('calc', i), ('exp', j)) or the reverse
        node1, node2 = pair
        if node1[0] == 'calc':
            calc_index = node1[1]
            exp_index = node2[1]
        else:
            calc_index = node2[1]
            exp_index = node1[1]

        matched_pair_indices.add((calc_index, exp_index))

    # Now, use these indices to access values in calc and exp
    matched_pairs = [(calc[i], exp[j]) for i, j in matched_pair_indices]
    assigned = {}
    calc_matched, exp_matched = [], []
    for calc_shift, exp_shift in matched_pairs:
        calc_matched.append(calc_shift)
        exp_matched.append(exp_shift)
        assigned[calc_shift] = exp_shift

    exclude_calc = list_difference(calc, calc_matched)
    exclude_exp = list_difference(exp, exp_matched)

    return assigned, exclude_exp, exclude_calc