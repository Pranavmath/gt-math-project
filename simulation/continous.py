import numpy as np
import json
import pprint
import scipy.integrate
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from tqdm import tqdm
import shutil 


# Define constants
UNLIMITED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]
ENZYME_CONCENTRATION = 100  # µM
INITIAL_METABOLITE_CONCENTRATION = 1500

TINKERED_REACTIONS = [
    {"kcat": 5.7, "keq": "+i", "products": {"FRUCTOSE-6P": 1, "Pi": 1}, "substrate_kms": {"FRUCTOSE-16-DIPHOSPHATE": 52.5}, "substrates": {"FRUCTOSE-16-DIPHOSPHATE": 1, "WATER": 1}},
    {"kcat": 1.0, "keq": "+i", "products": {"AMMONIA": 1, "GLT": 1, "PROTON": 1}, "substrate_kms": {"GLN": 1700.0}, "substrates": {"GLN": 1, "WATER": 1}},
]

with open("dataset.json", "r") as f:
    # List of reaction dictionaries.
    data = json.load(f)


# intial metabolite_concentrations = smth µM
metabolite_to_idx = {}

idx = 0

for rxn in data:
    for substrate in rxn["substrates"].keys():
        if substrate not in metabolite_to_idx.keys():
            metabolite_to_idx[substrate] = idx
            idx += 1
    for product in rxn["products"].keys():
        if product not in metabolite_to_idx.keys():
            metabolite_to_idx[product] = idx
            idx += 1


# no duplicate indexs
assert len(set(metabolite_to_idx.values())) == len(metabolite_to_idx.values())


idx_to_metabolite = {v: k for k,v in metabolite_to_idx.items()}


curr_metabolite_concentrations = np.full((len(metabolite_to_idx.values()), ), INITIAL_METABOLITE_CONCENTRATION)

#print(metabolite_to_idx)

last_t = 0
perturbed = 0

def calculate_concentration_changes(t, metabolite_concentrations, pbar, bruv, set_change):
    """
    Calculate the change in concentrations of metabolites based on reaction data.
    
    Args:
        metabolite_concentrations (np array): Initial concentrations of metabolites.
        
    Returns:
        dict: Changes in concentrations for each metabolite.
    """

    global last_t
    global perturbed

    # Initialize change in concentrations
    concentration_changes = np.zeros(metabolite_concentrations.shape)

    # Loop through each reaction
    for reaction in data:
        kcat = reaction["kcat"]
        keq = reaction["keq"]
        substrates = reaction["substrates"]
        products = reaction["products"]
        substrate_kms = reaction.get("substrate_kms", {})
        product_kms = reaction.get("product_kms", {})


        # Initialize reaction rates
        v_forward = 0
        v_reverse = 0

        vmax = ENZYME_CONCENTRATION * kcat


        # Calculate forward rate if substrate_kms exists
        if substrate_kms:
            vmax_forward = ENZYME_CONCENTRATION * kcat
            v_forward = vmax_forward
            for substrate, stoich in substrates.items():
                if substrate not in UNLIMITED_SUBSTRATES:
                    v_forward *= (metabolite_concentrations[metabolite_to_idx[substrate]] ** stoich) / (
                        substrate_kms[substrate] + metabolite_concentrations[metabolite_to_idx[substrate]]
                    )
            
            #v_forward = min(v_forward, vmax_forward)

        # Calculate reverse rate if product_kms exists
        if product_kms:
            vmax_reverse = ENZYME_CONCENTRATION * kcat
            v_reverse = vmax_reverse
            for product, stoich in products.items():
                if product not in UNLIMITED_SUBSTRATES:
                    v_reverse *= (metabolite_concentrations[metabolite_to_idx[product]] ** stoich) / (
                        product_kms[product] + metabolite_concentrations[metabolite_to_idx[product]]
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
                concentration_changes[metabolite_to_idx[substrate]] += change
            
            

        for product, stoich in products.items():
            if product not in UNLIMITED_SUBSTRATES:
                change = v_net * stoich
                concentration_changes[metabolite_to_idx[product]] += change
    
    pbar.update(round(t - last_t, 2))
    last_t = t

    
    if (2000 <= t) and (perturbed < 15000):
        #bruv = "DUTP"
        concentration_changes[metabolite_to_idx[bruv]] = set_change - metabolite_concentrations[metabolite_to_idx[bruv]]
        #print(concentration_changes[metabolite_to_idx["GMP"]])
        #print(t)
        perturbed += 1
    

    # make sure it doesn't go below 0
    for i in range(len(concentration_changes)):
        change = concentration_changes[i]
        conc = metabolite_concentrations[i]
        if (change < 0) and (conc + change < 0):
            concentration_changes[i] = -conc
    

    return concentration_changes


with open('simulation.npy', 'rb') as f:
    original = np.load(f)


pairs = [("Pi", 0), ("Pi", 1e+8), ("DUMP", 0), ("DUMP", 5000), ("S-ADENOSYLMETHIONINE", 0), ("S-ADENOSYLMETHIONINE", 5000)]


for idx, (bruv, set_change) in enumerate(pairs):
    NUM_ITER = 10000
    t = np.linspace(0, NUM_ITER, NUM_ITER * 10)


    # """
    with tqdm(total=NUM_ITER, unit="‰") as pbar:
        sol = scipy.integrate.solve_ivp(calculate_concentration_changes, [0, NUM_ITER], curr_metabolite_concentrations, method="LSODA", dense_output=True, args=[pbar, bruv, set_change])
    
    #print(sol)

    z = sol.sol(t)
    # """



    """
    with open("simulation.npy", "wb") as f:
        np.save(f, z)
    """




    """"
    datapoints = z[metabolite_to_idx["GMP"]]
    plt.plot(t, datapoints)
    plt.show()
    """

    print("finished sim")

    # """
    for metabolite in metabolite_to_idx.keys():
        datapoints = z[metabolite_to_idx[metabolite]]
        og = original[metabolite_to_idx[metabolite]]

        # Plot the iterative data with markers and a solid line
        plt.plot(t, og, label="Original Continous", color="purple", linestyle='-', linewidth=1.5)
        
        # Plot the continuous data with a dashed line and no markers
        plt.plot(t, datapoints, label="Perturbed Continous", color="orange", linestyle='-', linewidth=1.5)

        # Add title and axis labels with LaTeX for µM symbol
        plt.title(f"Comparison for {metabolite}", fontsize=16, fontweight='bold')
        plt.xlabel("Iteration", fontsize=14)
        plt.ylabel(r'Concentration ($\mu$M)', fontsize=14)

        # Add a legend with improved positioning and font size
        plt.legend(loc="upper left", fontsize=12, frameon=False)

        # Enable grid for better readability of the plot
        plt.grid(True, linestyle='--', alpha=0.7)

        # Improve the layout to avoid overlap
        plt.tight_layout()

        # Save the plot with higher resolution (300 dpi)
        plt.savefig(f"graphs/{metabolite}", dpi=600)
        
        # Clear figure for the next plot
        plt.clf()
    # """

    print("finished saving")

    shutil.make_archive(f"graphs{idx}", "zip", "graphs")

    print("finished zipping")

    last_t = 0
    perturbed = 0

