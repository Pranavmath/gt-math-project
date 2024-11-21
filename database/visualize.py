import json
import networkx as nx
from matplotlib import pyplot as plt

NOT_INCLUDED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]

G = nx.Graph()

with open("dataset.json", "r") as f:
    data = json.load(f)


for rxn in data:
    kcat = rxn["kcat"]
    substrate_kms = rxn.get("substrate_kms", {})
    product_kms = rxn.get("product_kms", {})


    for substrate in rxn["substrates"].keys():
        if substrate not in NOT_INCLUDED_SUBSTRATES:
            G.add_node(substrate)

            for product in rxn["products"].keys():
                G.add_node(product)
                    
                G.add_edge(substrate, product)


"""
with open("database/reactions.dat", "r") as f:
    temp = f.read()
    parsed_reactions = parse_reactions(temp)



for rxn in parsed_reactions.values():
    for substrate in rxn["left_substrates"].keys():
        G.add_node(substrate)

        for product in rxn["right_products"].keys():
            G.add_node(product)
                    
            G.add_edge(substrate, product)

"""


print("num of dsu: " + str(len(list(nx.connected_components(G)))))

nx.draw(G, with_labels=False, node_size=1)
plt.savefig("my_graph.pdf", dpi=1000)