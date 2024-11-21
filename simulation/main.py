import numpy as np
import json
import pprint
import seaborn as sns
import matplotlib.pyplot as plt


# Define constants
UNLIMITED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]
ENZYME_CONCENTRATION = 1  # µM
INITIAL_METABOLITE_CONCENTRATION = 1500

with open("dataset.json", "r") as f:
    # List of reaction dictionaries.
    data = json.load(f)


# intial metabolite_concentrations = smth µM
curr_metabolite_concentrations = {}

for rxn in data:
    for substrate in rxn["substrates"].keys():
        curr_metabolite_concentrations[substrate] = INITIAL_METABOLITE_CONCENTRATION
    for product in rxn["products"].keys():
        curr_metabolite_concentrations[product] = INITIAL_METABOLITE_CONCENTRATION



#print(len(data))

kcats = []
stoichs = []
kms = []


num_passes_filter = 0

def calculate_concentration_changes(metabolite_concentrations):
    """
    Calculate the change in concentrations of metabolites based on reaction data.
    
    Args:
        metabolite_concentrations (dict): Initial concentrations of metabolites.
        
    Returns:
        dict: Changes in concentrations for each metabolite.
    """

    global num_passes_filter

    # Initialize change in concentrations
    concentration_changes = {metabolite: 0 for metabolite in metabolite_concentrations}

    # Loop through each reaction
    for reaction in data:
        kcat = reaction["kcat"]
        keq = reaction["keq"]
        substrates = reaction["substrates"]
        products = reaction["products"]
        substrate_kms = reaction.get("substrate_kms", {})
        product_kms = reaction.get("product_kms", {})


        rxn_kms = []

        for km in substrate_kms.values():
            kms.append(km)
            rxn_kms.append(km)
        
        for km in product_kms.values():
            kms.append(km)
            rxn_kms.append(km)

        kcats.append(kcat)

        passes_filter = (substrate_kms or product_kms) and (max(rxn_kms) < 5000) and (kcat < 300)

        # -------------------------------------------------------

        if passes_filter:
            # Initialize reaction rates
            v_forward = 0
            v_reverse = 0

            vmax = ENZYME_CONCENTRATION * kcat

            num_passes_filter += 1

            # Calculate forward rate if substrate_kms exists
            if substrate_kms:
                vmax_forward = ENZYME_CONCENTRATION * kcat
                v_forward = vmax_forward
                for substrate, stoich in substrates.items():
                    if substrate not in UNLIMITED_SUBSTRATES:
                        v_forward *= (metabolite_concentrations[substrate] ** stoich) / (
                            substrate_kms[substrate] + metabolite_concentrations[substrate]
                        )
                
                #v_forward = min(v_forward, vmax_forward)

            # Calculate reverse rate if product_kms exists
            if product_kms:
                vmax_reverse = ENZYME_CONCENTRATION * kcat
                v_reverse = vmax_reverse
                for product, stoich in products.items():
                    if product not in UNLIMITED_SUBSTRATES:
                        v_reverse *= (metabolite_concentrations[product] ** stoich) / (
                            product_kms[product] + metabolite_concentrations[product]
                        )
                
                #v_reverse = min(v_reverse, vmax_reverse)

            # Handle equilibrium conditions
            if keq == "+i":
                v_net = v_forward  # Only forward reaction
            elif keq == "-i":
                v_net = -v_reverse  # Only reverse reaction
            else:
                # General equilibrium case
                keq_value = float(keq)
                if v_forward != 0:  # If forward rate is calculated, calculate reverse rate
                    v_reverse = v_forward / keq_value
                elif v_reverse != 0:  # If reverse rate is calculated, calculate forward rate
                    v_forward = v_reverse * keq_value

                v_net = v_forward - v_reverse

            # -vmax < v_net < vmax
            v_net = min(v_net, vmax)
            v_net = max(v_net, -vmax)

            # Update changes in metabolite concentrations (only for no unlimted ones)
            for substrate, stoich in substrates.items():
                if substrate not in UNLIMITED_SUBSTRATES:
                    change = -v_net * stoich
                    concentration_changes[substrate] += change
                
                stoichs.append(stoich)
                

            for product, stoich in products.items():
                if product not in UNLIMITED_SUBSTRATES:
                    change = v_net * stoich
                    concentration_changes[product] += change
                
                stoichs.append(stoich)
            
    return concentration_changes

NUM_ITER = 20000

datapoints = []

for _ in range(NUM_ITER):
    datapoints.append(curr_metabolite_concentrations["L-PHOSPHATIDATE"])

    for metabolite, change in calculate_concentration_changes(curr_metabolite_concentrations).items():
        curr_metabolite_concentrations[metabolite] += change

pprint.pprint(curr_metabolite_concentrations)

print(f"max {max(curr_metabolite_concentrations.values())}")
print(f"min {min(curr_metabolite_concentrations.values())}")

print(num_passes_filter)

plt.plot(datapoints)
plt.show()

"""
print(max(stoichs))
print(max(kcats))


asd = [i for i in kcats if i < 300]

print(len(asd)/len(kcats))


bsd = [i for i in kms if i < 5000]

print(len(bsd)/len(kms))


plt.hist(bsd, bins=100)
plt.show()
"""