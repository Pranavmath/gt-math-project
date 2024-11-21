import json
import networkx as nx
import pprint

NOT_INCLUDED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]

G = nx.Graph()

with open("rawdataset.json", "r") as f:
    data = json.load(f)

num = 0

filtered_data = []

for rxn in data:
    kcat = rxn["kcat"]
    substrate_kms = rxn.get("substrate_kms", {})
    product_kms = rxn.get("product_kms", {})

    rxn_kms = []

    for km in substrate_kms.values():
        rxn_kms.append(km)
    
    for km in product_kms.values():
        rxn_kms.append(km)
    
    # ----------------------------------------------


    passes_filter = (substrate_kms or product_kms) and (max(rxn_kms) < 5000) and (kcat < 300)

    if passes_filter:
        filtered_data.append(rxn)

        num += 1

        for substrate in rxn["substrates"].keys():
            if substrate not in NOT_INCLUDED_SUBSTRATES:
                G.add_node(substrate)

                for product in rxn["products"].keys():
                    G.add_node(product)
                        
                    G.add_edge(substrate, product)


print(num)

main_graph = max(nx.connected_components(G), key=lambda a: len(a))

processed_data = []

for rxn in filtered_data:
    substrate_kms = rxn.get("substrate_kms", {})
    product_kms = rxn.get("product_kms", {})
    substrates = rxn["substrates"]
    products = rxn["products"]

    assert substrate_kms or product_kms
    assert substrates and products

    in_main_graph = True

    for substrate in substrates:
        if substrate not in main_graph:
            in_main_graph = False
    for product in products:
        if product not in main_graph:
            in_main_graph = False
    
    
    if in_main_graph:
        processed_data.append(rxn)


pretty_json_str = pprint.pformat(processed_data, compact=True).replace("'",'"')

with open('dataset.json', 'w') as f:
    f.write(pretty_json_str)
