import json
import networkx as nx
from matplotlib import pyplot as plt
import plotly.graph_objects as go
import re
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

NOT_INCLUDED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]

G = nx.Graph()
concentrations = [3.66190913e-07, 1500.0, 3000.0, 312384.944, 6.4e-323, 21663.5982, 1500.0, 5843.52788, 1910.70933, 995.132038, 1.83e-322, 1.30017309e-143, 4231.24107, 899.960842, 1.01001369e-46, 29366.0796, 3000.0, 0, 1484.71941, 0, 6855.52029, 1515.28059, 4.01604528e-26, 504.656401, 276000.0, 0.0, 0.0, 0, 0, 0.0, 1500.0, 1597.28317, 0, 6000.0, 0, 3098.12028, 3.28433465e-36, 9916.54808, 16400.0117, 0, 3703.42036, 3.80568828, 2445.19784, 0.131106006, 1062.62073, 1500.0, 1937.37927, 1500.0, 70130.8266, 0.0, 0, 3151.01586, 0, 3000.0, 3000.0, 0.0, 1271500.0, 0.0, 0.0, 1.04861541e-52, 3583.12266, 2.5e-323, 2.22950269e-09, 1.0835504e-223, 3000.0, 3.5e-323, 1500.0, 997.183906, 0, 4031.22582, 1458.24222, 1636.8162, 1541.75778, 0, 0, 3000.0, 1500.0, 1500.0, 0.0, 0, 0, 0, 4500.0, 1408.24469, 1591.75531, 3000.0, 1166.56592, 0, 3.34838981e-156, 3000.0, 999.883296, 2999.65454, 0, 3125.17308, 1.47534281e-37, 9.24483139, 0, 3000.0, 6.23056647e-05, 6.23056647e-05, 3814.75964, 3000.40805, 6.78084154e-07, 1500.0, 1500.0, 5.07657349e-177, 1.69955459e-84, 3000.0, 0, 3000.0, 4118.6771, 0.23526799, 2999.76473, 2.08748807e-50, 3000.0, 4.83785149e-50, 3000.0, 0, 3729.9184, 0, 510.238646, 4718.46823, 11090.1173, 1677.48894, 70280.4344, 2151500.0, 1456.97716, 1543.02284, 0, 0, 1500.0, 3000.0, 0, 3000.0, 5.6611615e-72, 4506.09794, 93.3980487, 5.12700511, 1583.0767, 2.76394138e-18, 3000.0, 1487.11156, 1512.88844, 5.92354404e-35, 1.03995362e-112, 3000.0, 1500.0, 1500.0, 3.6085248e-261, 1636.64744, 0.193502752, 4502.01404, 1464.08469, 1535.91531, 1.35839387e-19, 3000.0, 8416.54817, 0, 4568.38615, 0.0, 0, 3000.0, 0, 3000.0, 3156.47212, 4343.52788, 0, 3000.0, 0, 2510.77259, 2e-322, 1219.3548, 1780.6452, 0, 4545.11021, 1500.0, 1078.08149, 0.115799004, 6.02352981e-271, 1500.0, 1.52506394e-107, 0, 3049.03479, 0, 3000.0, 1976.83376, 1523.28294, 1523.28294, 999.883296, 0, 0.0, 2149.61177, 1499.99991, 3000.0, 3000.0, 1.05261838e-07, 3012.28916, 1.1e-322, 0, 1.53697155e-09, 405.83998, 0.881939715, 0.881939715, 1.43026792e-10, 3000.0, 4.77992096e-31, 3000.0, 1499.50017, 1500.49983, 4100.0, 4.50977339e-97, 0, 4500.0, 15.7454815, 1.69384473e-211, 3000.0, 0.345462289, 2999.65454, 0, 0.0, 2824.62027, 5e-324, 3000.0, 7.10408244e-212, 1.3e-322, 2355.03437, 7.25235788e-198, 5996.06321, 5.27038806e-241, 1995.17638, 1484.81495, 1515.18505, 0.0, 56500.0, 3.68401116e-112, 0, 0, 3074.29632, 7.53842948e-158, 4749.72237, 2.02562903e-226, 3000.0, 3.64825907e-19, 3000.0, 3000.0, 0, 3014.01063, 0, 3000.0, 5.28628656e-272, 3850.38823, 3101.59405, 1442.94285, 1557.05715, 9.40989696e-05, 9.40989696e-05, 1500.0, 113.685058, 2886.31494, 2100.03916, 2100.03916, 0, 3000.0, 5.61829643e-43, 3000.0, 681500.0, 2.48645021e-187, 0.0164785017, 2249.99176, 9.62656903e-16, 6000.0, 1.7648e-318, 1500.0, 194.812287, 2805.18771, 1397.34594, 1602.65406, 191500.0, 771.293124, 2898.40595, 1500.0, 5.34672318e-279, 2762.86378, 8.74815285e-07, 3000.0, 0, 3000.0, 3.1e-322, 0, 1500.0, 0.000961948787, 0.000961948787, 2999.99904, 1.70949628e-176, 3000.0, 0.0, 0, 3000.0, 2863.35256, 8.30125491e-152, 0, 3002.39502, 1.79871886e-18, 4093.27808, 0, 993.194371, 1500.0, 1500.0, 1500.0, 1499.56461, 1500.43539, 3053.75108, 0, 2159.76029, 840.239705, 0.0, 141500.0, 1500.0, 1.24e-322, 5955.28325, 0, 3000.0, 1500.0, 1611.96314, 0, 44.7167515, 1.10955472e-251, 3000.0, 3000.0, 3.26537106e-97, 3000.0, 1544.50577, 3164.96433, 1500.0, 1500.0, 280.645197, 2719.3548, 5.60292566e-15, 4401.47495, 0, 3000.0, 685.648402, 685.648402, 2314.3516, 3.01035903e-07, 3000.0, 12500.0, 41027.5, 0.0, 0, 3000.0, 5.12210516e-124, 3000.0, 0, 0, 2.75294283e-16, 3000.0, 3.5802577e-119]

epsilon = 1e-9  # A small value to avoid log(0)
log_concentrations = np.log10(np.array(concentrations) + epsilon)

# Normalize the log-transformed concentrations to range [0, 1]
norm = Normalize(vmin=log_concentrations.min(), vmax=log_concentrations.max())

# Use a perceptually uniform colormap (e.g., 'coolwarm')
cmap = plt.cm.coolwarm

# Map the normalized concentration values to colors
sm = ScalarMappable(cmap=cmap, norm=norm)
node_colors = sm.to_rgba(log_concentrations, bytes=True)

# Convert to rgba format for Plotly
node_colors = [f'rgba({r},{g},{b},{a/255})' for r, g, b, a in node_colors]


def parse_reactions(data):
    # Split reactions by `//`
    reactions = data.strip().split("//")
    result = {}

    for reaction in reactions:
        if not reaction.strip():
            continue
        
        # Extract UNIQUE-ID
        unique_id_match = re.search(r"UNIQUE-ID - ([^\n]+)", reaction)
        unique_id = unique_id_match.group(1).strip() if unique_id_match else None

        # Extract GIBBS-0
        gibbs_match = re.search(r"GIBBS-0 - ([^\n]+)", reaction)
        gibbs_0 = float(gibbs_match.group(1).strip()) if gibbs_match else None

        # extract ec number
        ec_match = re.search(r"EC-NUMBER - ([^\n]+)", reaction)
        ec_number = ec_match.group(1).strip() if ec_match else None

        direction_hash = {
            "PHYSIOL-LEFT-TO-RIGHT": 0,
            "LEFT-TO-RIGHT": 0,
            "PHYSIOL-RIGHT-TO-LEFT": 1,
            "RIGHT-TO-LEFT": 1,
            "REVERSIBLE": 2
        }

        # Extract REACTION-DIRECTION
        direction_match = re.search(f"REACTION-DIRECTION - ({'|'.join(direction_hash.keys())})", reaction)

        direction = direction_hash[direction_match.group(1)] if direction_match else None

        # Extract LEFT and RIGHT substances
        left_matches = re.findall(r"LEFT - ([^\n]+)", reaction)
        right_matches = re.findall(r"RIGHT - ([^\n]+)", reaction)

        # Extract coefficients for LEFT and RIGHT substances
        left_coefficients = {}
        for match in re.finditer(r"LEFT - ([^\n]+)(\n\^COEFFICIENT - (\d+))?", reaction):
            substance = match.group(1).strip()
            coefficient = int(match.group(3)) if match.group(3) else 1
            left_coefficients[substance] = coefficient

        right_coefficients = {}
        for match in re.finditer(r"RIGHT - ([^\n]+)(\n\^COEFFICIENT - (\d+))?", reaction):
            substance = match.group(1).strip()
            coefficient = int(match.group(3)) if match.group(3) else 1
            right_coefficients[substance] = coefficient

        # add parsed reaction data
        result[unique_id] = {
            "unique_id": unique_id,
            "gibbs_0": gibbs_0,
            "direction": direction,
            "left_substrates": left_coefficients,
            "right_products": right_coefficients,
            "original_text": reaction,
            "ec_number": ec_number
        }
    
    return result


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

coloring = nx.coloring.greedy_color(G, strategy="largest_first")
chromatic_number = len(set(coloring.values()))
print(chromatic_number)

print(len(G.nodes))
print(len(G.edges))

"""
node_colors = ["skyblue"] * len(G.nodes())  # Default color for all nodes
rss_nodes = [
    "CPD-15384",
    "CPD-20750",
    "L-RIBULOSE-5-P",
    "CPD-7248",
    "CPD-7249",
    "CPD-381",
    "HYDROGEN-PEROXIDE",
    "LACTALD",
    "LL-DIAMINOPIMELATE",
    "R-2-HYDROXYGLUTARATE",
    "ADP",
    "2-PG",
    "D-GLT",
    "GALACTOSE-1P",
    "CPD0-1445",
    "CPD0-1965",
    "ETOH",
    "CARBON-DIOXIDE",
    "NADH",
    "NAD",
    "PROTON",
    "WATER",
    "DEHYDRO-DEOXY-GALACTONATE-PHOSPHATE",
    "CPD-20746"
]
iss_nodes = [
    "DIHYDRONEOPTERIN-P", "PYRAZINOIC-ACID", "ETHANOL-AMINE", "2-DEHYDRO-3-DEOXY-D-GALACTONATE",
    "ACETOL", "7-METHYLGUANOSINE-5-PHOSPHATE", "UDP-ACETYL-CARBOXYVINYL-GLUCOSAMINE", "4-CYTIDINE-5-DIPHOSPHO-2-C",
    "1-AMINO-PROPAN-2-ONE-3-PHOSPHATE", "CPD-25028", "CPD-16720", "SO3", "OXIDIZED-GLUTATHIONE", "GLYCOLLATE", 
    "FADH2", "OCTAPRENYL-DIPHOSPHATE", "16S-rRNA-N3-methyluracil1498", "CPD-389", "D-ALPHABETA-D-HEPTOSE-7-PHOSPHATE",
    "MALTOHEXAOSE", "CPD-560", "6-Dimethylallyladenosine37-tRNAs", "2-3-DIHYDROXYBENZOATE", "MANNOSE-1P", 
    "RIBOFLAVIN", "CPD-12587", "ADENOSINE_DIPHOSPHATE_RIBOSE", "ADENYLOSUCC", "PYRIDOXAL", "PORPHOBILINOGEN", 
    "GLY", "2-DEOXY-D-GLUCOSE", "3-HYDROXYBENZOATE", "KDO-8P", "CL-", "B-ALANINE", "3OH-4P-OH-ALPHA-KETOBUTYRATE", 
    "DIHYDROXYNAPHTHOATE", "CPD0-2480", "S-24-DINITROPHENYLGLUTATHIONE", "FUM", "DNA-Containing-N6-Methyladenine", 
    "DTDP-DEOH-DEOXY-GLUCOSE", "XANTHINE", "4-AMINO-4-DEOXYCHORISMATE", "4-AMINO-BUTYRATE", "CPD-782", 
    "S-ADENOSYLMETHIONINAMINE", "ADENOSINE", "UDP-MURNAC-TETRAPEPTIDE", "DEHYDRO-3-DEOXY-L-RHAMNONATE", "UMP", 
    "BETA-D-FRUCTOSE", "23S-rRNA-2-O-methyluridine2552", "DIAMINONONANOATE", "S-ADENOSYL-4-METHYLTHIO-2-OXOBUTANOATE", 
    "8-AMINO-7-OXONONANOATE", "Carboxylates", "CPD0-932", "5-P-RIBOSYL-N-FORMYLGLYCINEAMIDE", 
    "3-DEOXY-D-ARABINO-HEPTULOSONATE-7-P", "UDP", "MAL", "D-BETA-D-HEPTOSE-1-P", "PROPIONYL-COA", "DEOXYGUANOSINE", 
    "*CPD-14447", "GLYCEROL-3P", "DIHYDROFOLATE-GLU-N", "CMP-KDO", "SORBITOL", "DIAMINO-OH-PHOSPHORIBOSYLAMINO-PYR", 
    "PHENYLACETATE", "MYO-INOSITOL", "GDP-4-DEHYDRO-6-DEOXY-D-MANNOSE", "HS", "GLYCEROL", "CPD-13243", 
    "D-MYO-INOSITOL-1-MONOPHOSPHATE", "23S-rRNA-N3-methylpseudouridine1915", "tRNA-Containing-N1-Methylguanine-37", 
    "DEOXYURIDINE", "INDOLE", "DIHYDROXYPENTANEDIONE", "PANTETHEINE-P", "G3P", "HYPOXANTHINE", "CIT", "CPD-14443", 
    "AMP", "L-DI-GMP", "Ox-Thioredoxin", "ADENOSYL-HOMO-CYS", "CPD0-2190", "CPD-9646", "ACET", "D-ALANINE", 
    "tRNA-containing-5Me-uridine54", "4-hydroxybenzoate", "CARBAMATE", "L-PHOSPHATIDATE", "2K-4CH3-PENTANOATE", 
    "D-Ribofuranose", "CADAVERINE", "AMMONIA", "CPD-9923", "DEOXYXYLULOSE-5P", "FORMYL-COA", "tRNA-Arg-inosine34", 
    "TMP", "AMMONIUM", "MET", "tRNAs-containing-epoxy-quenosine", "D-LACTATE", "DEOXYCYTIDINE"
]
dli_nodes = [
    "Dodecanoyl-ACPs", "MALONYL-ACP", "CPD-13211", "CPD-558", "CPD-19179", "L-ALPHA-ALANINE", "11-DEOXYCORTICOSTERONE"
]
ili_nodes = [
    "CPD-16569", "GLUCOSAMINE-1P", "3-oxo-myristoyl-ACPs", "DEOXYADENOSINE", "3-ENOLPYRUVYL-SHIKIMATE-5P", 
    "CPD-618", "P-NITROPHENOL", "DEOXY-RIBOSE-5P", "CYS-GLY", "5-OXOPROLINE", "GLC-6-P", "THYMINE", "TARTRATE", 
    "2-DH-3-DO-D-ARABINONATE", "D-XYLONATE", "IMP", "MANNITOL-1P", "PRECURSOR-Z", "CARBAMYUL-L-ASPARTATE", 
    "CPD-16843", "UDP-N-ACETYL-D-GLUCOSAMINE", "5-KETO-4-DEOXY-D-GLUCARATE"
]
di_nodes = [
    "DAMP"
]
irr_nodes = [
    "XANTHOSINE-5-PHOSPHATE", "SER", "ERYTHRONATE-4P", "Glucopyranose", "P3I", "CPD0-930", "D-SEDOHEPTULOSE-7-P", 
    "SUC", "METHYLENE-THF-GLU-N", "GUANOSINE-5DP-3DP", "CYTIDINE", "ACETYL-COA", "OXALACETIC_ACID", "METHYL-GLYOXAL", 
    "HOMO-CYS", "D-GLUCOSAMINE-6-P", "THF-GLU-N", "L-ASPARTATE", "PYRUVATE", "DUMP", "GAP", "LYS", "FORMYL-THF-GLU-N", 
    "GMP", "GLUTATHIONE"
]
dss_nodes = [
    "ENOL-OXALOACETATE", "GLYCOLALDEHYDE", "PYRAZINAMIDE", "L-1-GLYCEROPHOSPHORYLETHANOL-AMINE", "3-P-SERINE", 
    "DIHYDRONEOPTERIN-P3", "PHENYLACETALDEHYDE", "4-PHOSPHONOOXY-THREONINE", "D-GALACTONATE", "CPD0-2030", 
    "2-C-METHYL-D-ERYTHRITOL-4-PHOSPHATE", "FAD", "PAPS", "CPD0-2511", "CPD-67", "SHIKIMATE-5P", "CPD-8887", 
    "CPD0-1202", "CPD0-2552", "2-KETO-3-DEOXY-6-P-GLUCONATE", "16S-rRNA-uracil1498", "CPD-8989", "CPD-1091", 
    "CPD-397", "D-SERINE", "CPD-4211", "D-CYSTEINE", "DIHYDRO-DIOH-BENZOATE", "UTP", "CPD0-1133", "tRNA-adenine-37", 
    "ACETOACETYL-COA", "5-METHYLTHIOADENOSINE", "CPD-14762", "5-AMINO-LEVULINATE", "HCO3", "L-CYSTATHIONINE", 
    "3-P-HYDROXYPYRUVATE", "DI-H-URACIL", "FMN", "PYRIDOXAL_PHOSPHATE", "CPD0-929", "2-DEOXY-D-GLUCOSE-6-PHOSPHATE", 
    "CPD-9925", "ADP-L-GLYCERO-D-MANNO-HEPTOSE", "CPD0-2479", "1-CHLORO-24-DINITROBENZENE", "CPD-264", "CPD-18118", 
    "CPD0-2184", "34-DIHYDROXYPHENYLACETYL-COA", "DNA-Adenines", "CPD0-2461", "CPD-660", "CPD-448", "DTDP-D-GLUCOSE", 
    "CPD-606", "DIHYDRO-THYMINE", "L-ALA-GAMMA-D-GLU-DAP", "4-AMINO-BUTYRALDEHYDE", "FRU1P", "FRUCTOSE-16-DIPHOSPHATE", 
    "L-RHAMNONATE", "CPD-3706", "CPD-66", "XYLULOSE-5-PHOSPHATE", "23S-rRNA-uridine-2552", "TAGATOSE-1-6-DIPHOSPHATE", 
    "C1", "DIHYDROXY-ACETONE-PHOSPHATE", "RIBOSE-5P", "OROTIDINE-5-PHOSPHATE", "Acyl-Phosphates", "DGTP", "GLYOX", 
    "23-Diaminopropanoate", "ERYTHROSE-4P", "PHOSPHO-ENOL-PYRUVATE", "5-PHOSPHO-RIBOSYL-GLYCINEAMIDE", 
    "MESO-DIAMINOPIMELATE", "KDO", "2-KETOGLUTARATE", "D-SORBITOL-6-P", "DELTA3-ISOPENTENYL-PP", "FARNESYL-PP", 
    "CPD-157", "METHYL-MALONYL-COA", "GDP-TP", "D-BETA-D-HEPTOSE-17-DIPHOSPHATE", "DGMP", "CTP", "RIBOSE-1P", 
    "GLYCEROPHOSPHOGLYCEROL", "GERANYL-PP", "23S-rRNA-pseudouridine1915", "CYS", "DATP", "GDP-MANNOSE", 
    "R-4-PHOSPHOPANTOTHENOYL-L-CYSTEINE", "CPD-564", "L-ASPARTATE-SEMIALDEHYDE", "Guanine37-in-tRNA", "CPD-207", 
    "CPD-3711", "TRP", "3-5-ADP", "GLC-1-P", "Uracil-54-in-tRNA", "ISOCHORISMATE", "ACETALD", "C-DI-GMP", "D-GALACTARATE", 
    "UNDECAPRENYL-DIPHOSPHATE", "N-ACETYL-D-GLUCOSAMINE-6-P", "Red-Thioredoxin", "CPD-8990", "CYTOSINE", "GTP", 
    "DIACYLGLYCEROL-PYROPHOSPHATE", "UDP-MANNAC", "FRUCTOSE-6P", "CPD0-2339", "CPD-548", "2-D-THREO-HYDROXY-3-CARBOXY-ISOCAPROATE",
    "DI-H-OROTATE", "ASN", "5-METHYL-THF-GLU-N", "CPD0-2331", "CHORISMATE", "OXALYL-COA", "DUTP", "CPD-26279", "GLT", 
    "CPD-9924", "URIDINE", "GLN", "ADENINE", "S-ADENOSYLMETHIONINE", "5-10-METHENYL-THF-GLU-N", "tRNA-with-7-aminomethyl-7-deazaguanine",
    "GUANINE", "S-LACTOYL-GLUTATHIONE", "TTP", "tRNA-Arg-adenosine34", "PRPP", "DCMP"
]
ii_nodes = [
    "CO-A", "FORMATE", "URACIL", "Pi", "PPI"
]

for node in rss_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "red"
for node in iss_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "purple"
for node in dli_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "blue"        
for node in ili_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "green"  
for node in irr_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "yellow"  
for node in dss_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "orange"  
for node in ii_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "pink"  
for node in di_nodes:
    if node in G.nodes():
        node_index = list(G.nodes()).index(node)  # Find the index of the node
        node_colors[node_index] = "brown"  
"""
pos = nx.spring_layout(G, dim=3, seed=42, k=1.5, iterations=200)
x_nodes = [pos[node][0] for node in G.nodes()]  # X-coordinates
y_nodes = [pos[node][1] for node in G.nodes()]  # Y-coordinates
z_nodes = [pos[node][2] for node in G.nodes()] 
edge_x = []
edge_y = []
edge_z = []
for edge in G.edges():
    x0, y0, z0 = pos[edge[0]]
    x1, y1, z1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])  # None to break the line
    edge_y.extend([y0, y1, None])
    edge_z.extend([z0, z1, None])
node_trace = go.Scatter3d(
    x=x_nodes,
    y=y_nodes,
    z=z_nodes,
    mode='markers',
    marker=dict(size=3.5, color=node_colors),  # Customize node size and color
    text=list(G.nodes()),  # Node labels
    hoverinfo='text'
)

# Create 3D scatter plot for edges
edge_trace = go.Scatter3d(
    x=edge_x,
    y=edge_y,
    z=edge_z,
    mode='lines',
    line=dict(width=0.5, color='black'),  # Edge style
    hoverinfo='none'
)

# Combine traces
fig = go.Figure(data=[edge_trace, node_trace])
fig.update_layout(
    showlegend=False,
    scene=dict(
        xaxis=dict(showbackground=False),
        yaxis=dict(showbackground=False),
        zaxis=dict(showbackground=False)
    ),
    title="Partial E. coli Metabolic Network"
)

fig.show()
plt.figure(figsize=(12, 12))  # Adjust the figure size for better spacing




print(len(data))

print("num of dsu: " + str(len(list(nx.connected_components(G)))))

nx.draw(G, with_labels=False, node_size=1)
plt.savefig("my_graph.pdf", dpi=1000)