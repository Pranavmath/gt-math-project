import re
from collections import defaultdict
import statistics
import json
import pprint

NOT_INCLUDED_SUBSTRATES = ["ADP", "ATP", "PROTON", "WATER", "NAD", "NADH", "NADP", "NAD+", "CO2", "OXYGEN-MOLECULE", "NADPH", "NADH-P-OR-NOP", "CARBON-DIOXIDE", "NADc-P-OR-NOP"]


def parse_enzrxns(file_content):
    # Split the file content by "//" to separate each enzyme reaction
    reactions = file_content.strip().split("//")
    
    # Prepare a list to hold parsed data
    parsed_reactions = []
    
    # Regular expressions for the required fields
    unique_id_pattern = re.compile(r"UNIQUE-ID - (\S+)")
    reaction_pattern = re.compile(r"REACTION - (\S+)")
    enzyme_pattern = re.compile(r"ENZYME - (\S+)")
    kcat_pattern = re.compile(r"KCAT - ([\d\.\-eE]+)")
    km_pattern = re.compile(r"KM - ([\d\.\-eE]+)\s*\^SUBSTRATE - (\S+)")

    for reaction in reactions:
        if not reaction.strip():
            continue
        
        # Extract data
        unique_id = unique_id_pattern.search(reaction)
        reaction_name = reaction_pattern.search(reaction)
        enzyme = enzyme_pattern.search(reaction)
        
        # Extract KCAT values (no substrate association needed)
        kcat_matches = kcat_pattern.findall(reaction)
        kcat_values = [float(kcat) for kcat in kcat_matches]
        
        # Extract KM values and their substrates
        km_matches = km_pattern.findall(reaction)
        km_values = defaultdict(list)
        for km, substrate in km_matches:
            # strip bar from substrate name

            substrate = substrate.strip(" ")
            substrate = substrate.strip("|")
            substrate = substrate.strip(" ")


            km_values[substrate].append(float(km))
        
        # Append parsed reaction data to the list
        parsed_reactions.append({
            "unique_id": unique_id.group(1) if unique_id else None,
            "reaction": reaction_name.group(1) if reaction_name else None,
            "enzyme": enzyme.group(1) if enzyme else None,
            "kcat_values": kcat_values,
            "substrate_kms": dict(km_values),
        })
    
    return parsed_reactions


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



with open("database/enzrxns.dat", "r") as file:
    data = file.read()
    enzrxns = parse_enzrxns(data)


with open("database/reactions.dat", "r") as f:
    data = f.read()
    parsed_reactions = parse_reactions(data)


with open("database/grace.json", "r") as f:
    grace = json.load(f)

def hasallKms(subsrateKms, needsTheseSubstrates):
    """
    subsrateKms - dictionary of each substrate and their km values
    needsTheseSubstrates - list of substrates we need
    """

    hasall = True

    for substrate in needsTheseSubstrates:
        if substrate not in NOT_INCLUDED_SUBSTRATES:
                # if the substrate is in the dictionary and has at least 1 km value for that substrate in the dictionary
                if (substrate in subsrateKms.keys()) and (len(subsrateKms[substrate]) > 0):
                    hasall = hasall and True
                else:
                    hasall = False
    
    return hasall


def getallKms(subsrateKms, needsTheseSubstrates):
    """
    subsrateKms - dictionary of each substrate and their km values
    needsTheseSubstrates - list of substrates we need
    hasallKms is true
    """

    assert hasallKms(subsrateKms, needsTheseSubstrates) == True

    km_substrates = {}

    for substrate in needsTheseSubstrates:
        if substrate not in NOT_INCLUDED_SUBSTRATES:
                km_substrates[substrate] = statistics.median(subsrateKms[substrate])
    
    return km_substrates


i = 0

all_data_points = []


for enzrxn in enzrxns:
    # only 1 null reaction now in parsed reactions lol
    if enzrxn["reaction"] in parsed_reactions.keys():
        parsed_reaction = parsed_reactions[enzrxn["reaction"]]

        """
        direction of -1 means we dont have all kms for products or reactants
        direction of 0 means we have it for all the substrates => forward
        direcrtion of 1 means we have it for all the products => reverse
        """
        direction = -1

        if hasallKms(enzrxn["substrate_kms"], list(parsed_reaction["left_substrates"].keys())):
            direction = 0
        if hasallKms(enzrxn["substrate_kms"], list(parsed_reaction["right_products"].keys())):
            direction = 1
        


        # no (parsed_reaction["gibbs_0"] != None) since we dont need it (grace yay!)
        good = (len(enzrxn["kcat_values"]) > 0) and (parsed_reaction["direction"] != None) and (parsed_reaction["gibbs_0"] != None)

        matching_direction = (direction == 0 and parsed_reaction["direction"] == 0) or (direction == 1 and parsed_reaction["direction"] == 1)
        

        """
        if (direction != -1) and good and (parsed_reaction["direction"] == 2):
            #print(parsed_reaction["gibbs_0"])
            print(parsed_reaction["ec_number"])

            reaction_str = f'{" + ".join(parsed_reaction["left_substrates"].keys())} <=> {" + ".join(parsed_reaction["right_products"].keys())}'

            print(reaction_str)
            print("----------------------------------------")
            i += 1
        """


        data = {}

        if good and matching_direction:
            if (direction == 0):
                data["substrate_kms"] = getallKms(enzrxn["substrate_kms"], list(parsed_reaction["left_substrates"].keys()))
                data["keq"] = "+i"
            
            if (direction == 1):
                data["product_kms"] = getallKms(enzrxn["substrate_kms"], list(parsed_reaction["right_products"].keys()))
                data["keq"] = "-i"
            
            data["kcat"] = statistics.median(enzrxn["kcat_values"])

            data["substrates"] = parsed_reaction["left_substrates"]
            data["products"] = parsed_reaction["right_products"]

            i += 1
            
        
        # currently before grace does her thing we are just assuming k_eq is 1
        if (parsed_reaction["direction"] == 2) and good and (direction != -1):
            if (direction == 0):
                data["substrate_kms"] = getallKms(enzrxn["substrate_kms"], list(parsed_reaction["left_substrates"].keys()))
            
            if (direction == 1):
                data["product_kms"] = getallKms(enzrxn["substrate_kms"], list(parsed_reaction["right_products"].keys()))

            # --------------------------

            keq = -1

            ec = parsed_reaction["ec_number"]

            reaction_str = f'{" + ".join(parsed_reaction["left_substrates"].keys())} <=> {" + ".join(parsed_reaction["right_products"].keys())}'


            if (ec in grace.keys()):
                keq = float(grace[ec])
            if (reaction_str in grace.keys()):
                keq = float(grace[reaction_str])
            
            assert keq != -1

            # ----------------------------

            data["keq"] = keq

            data["kcat"] = statistics.median(enzrxn["kcat_values"])

            data["substrates"] = parsed_reaction["left_substrates"]
            data["products"] = parsed_reaction["right_products"]

            i += 1


        # not empty
        if data:
            all_data_points.append(data)
    



pretty_json_str = pprint.pformat(all_data_points, compact=True).replace("'",'"')

with open('rawdataset.json', 'w') as f:
    f.write(pretty_json_str)




print(i)
print(len(enzrxns))
