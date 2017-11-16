import requests
import math

def get_ftrd():
    url = "http://coleykursarlab.chpc.utah.edu/api/featuretablerawdata/"
    ftrd = []
    count = 1

    while True:
        if count % 10 == 0:
            print(count, "pages processed")
            
        resp = requests.get(url)
        data = resp.json()
        ftrd += data["results"]
        
        if data["next"] is None:
            print(count, "pages total")
            break
        else:
            url = data["next"]

        count += 1
    
    return ftrd

# match features across samples from ftrd
def build_feature_table(raw_data, mz_error, rt_error):
    feature_table = {}
    current_feature_number = 0
    for datum in raw_data:
        pcid = datum["PC_ID"]
        mz = datum["MZ"]
        rt = datum["RT"]
        sample = datum["sample"]
        status = "not matched"
        for feature in feature_table.keys():
            if abs(feature_table[feature]["avg_mz"] - mz) <= mz_error and \
                abs(feature_table[feature]["avg_rt"] - rt) <= rt_error:
                feature_table[feature]["mz"] += [mz]
                feature_table[feature]["rt"] += [rt]
                feature_table[feature]["PC_IDs"] += [pcid]
                feature_table[feature]["sample"] += [sample]
                feature_table[feature]["avg_mz"] = sum(feature_table[feature]["mz"])/len(feature_table[feature]["mz"])
                feature_table[feature]["avg_rt"] = sum(feature_table[feature]["rt"])/len(feature_table[feature]["rt"])
                feature_table[feature]["TIC"] += [datum["TIC"]]
                status = "matched"
                #if pcid in feature_table[feature]:
                #    feature_table[feature][pcid] += datum["TIC"]
                #else:
                #    feature_table[feature][pcid] = datum["TIC"]
                break
        if status == "not matched":
            current_feature_number += 1
            feature_table[current_feature_number] = {
                "mz": [mz],
                "rt": [rt],
                "PC_IDs": [pcid],
                "sample": [sample],
                "TIC": [datum["TIC"]],
                #pcid: datum["TIC"],
                "avg_mz": mz,
                "avg_rt": rt
                }
    
    return feature_table


# build pc_id table 
def build_pc_id(feature_table):
    pc_id_table = {}
    
    for feature in feature_table.keys():
        for pcid in list(set(feature_table[feature]["PC_IDs"])):
            tics = [feature_table[feature]["TIC"][i] for i, x in enumerate(feature_table[feature]["PC_IDs"]) if x == pcid]
            if pcid in pc_id_table:
                pc_id_table[pcid]["features"] += [feature]
                pc_id_table[pcid]["TICs"] += [sum(tics) / len(tics)]
            else:
                pc_id_table[pcid] = {
                    "features": [feature],
                    "TICs": [sum(tics) / len(tics)]
                }
    return pc_id_table

# create function for cosine score
def cosine_score(compound_1, compound_2, abundance, features):
    tics_1 = []
    tics_2 = []
    for feature in compound_1['features']:
        if feature in compound_2['features']:
            tics_1 += [float(compound_1[abundance][compound_1[features].index(feature)])]
            tics_2 += [float(compound_2[abundance][compound_2[features].index(feature)])]
    if len(tics_1) > 0:
        cosine_score = sum([tics_1[x] * tics_2[x] for x in range(len(tics_1))]) / \
        ((sum([x ** 2 for x in compound_1[abundance]]) ** (0.5)) * (sum([x ** 2 for x in compound_2[abundance]]) ** (0.5)))
    else: 
        cosine_score = 0
    return cosine_score

# build compound table. Combine PC_IDs that match with a high enough cosine score into a compound.
def build_compound_table(pc_id_table, min_cos_score):
    compound_table_temp = {}
    sorted_pc_ids = sorted(pc_id_table.keys(), key = lambda k: len(pc_id_table[k]["features"]), reverse = True)
    current_compound_number = 0

    for pcid in sorted_pc_ids:
        status = "not_matched"
        for compound in compound_table_temp:
            if cosine_score(pc_id_table[pcid], compound_table_temp[compound], abundance="TICs", features="features") > min_cos_score:
                for feature in pc_id_table[pcid]["features"]:
                    if feature in compound_table_temp[compound]["features"]:
                        compound_table_temp[compound]["TICs"][compound_table_temp[compound]["features"].index(feature)] += \
                        pc_id_table[pcid]["TICs"][pc_id_table[pcid]["features"].index(feature)]
                    else:
                        compound_table_temp[compound]["features"] += [feature]
                        compound_table_temp[compound]["TICs"] += [pc_id_table[pcid]["TICs"][pc_id_table[pcid]["features"].index(feature)]]
                status = "matched"
                break
        if status == "not_matched":
            current_compound_number += 1
            compound_table_temp[current_compound_number] = {
                "features": pc_id_table[pcid]["features"],
                "TICs": pc_id_table[pcid]["TICs"]
            }
    
    return compound_table_temp

# build function to match filled features to compounds:
def fill_compounds(filled_features, compound_table, min_cos_score = 0.3):
    filled_compounds = {}

    for sample in filled_features:
        for compound in compound_table:
            shared_features = [x for x in filled_features[sample]["feature_number"] if x in compound_table[compound]["features"]]
            if float(len(shared_features)) / float(len(compound_table[compound]["features"])) > 0.1:
                sample_compound = {
                    "features": shared_features,
                    "TICs": [filled_features[sample]["TIC"][filled_features[sample]["feature_number"].index(x)] for x in shared_features]
                }
                if cosine_score(sample_compound, compound_table[compound], abundance = "TICs", features = "features") > min_cos_score:
                    if sample in filled_compounds:
                        filled_compounds[sample]["compound"] += [compound]
                        filled_compounds[sample]["TIC"] += [sum(sample_compound["TICs"])]
                    else:
                        filled_compounds[sample] = {
                            "compound": [compound],
                            "TIC": [sum(sample_compound["TICs"])]
                        }
                    filled_features[sample]["TIC"] = [x for i,x in enumerate(filled_features[sample]["TIC"]) if i not in \
                    [filled_features[sample]["feature_number"].index(k) for k in shared_features]]
                    filled_features[sample]["feature_number"] = [x for i,x in enumerate(filled_features[sample]["feature_number"]) if i not in \
                    [filled_features[sample]["feature_number"].index(k) for k in shared_features]]
    
    return filled_compounds
    

##### The rest of this is old code that may not be necessary anymore.
#for feature in pc_id:
#    current_pcid = feature["PC_ID"]
#    if current_pcid in temp4:
#        if feature["feature"] in temp4[current_pcid]["features"]:
#            temp4[current_pcid]["average_TIC"][temp4[current_pcid]["features"].index(feature["feature"])] += \
#                feature["average_TIC"]
#        else:
#            temp4[current_pcid]["features"] += [feature["feature"]]
#            temp4[current_pcid]["average_TIC"] += [feature["average_TIC"]]
#    else:
#        temp4[current_pcid] = {
#            "features": [feature["feature"]],
#            "average_TIC": [feature["average_TIC"]],
#            "compound_number": feature["compound_number"]
#        }
#
#for k in sorted(temp4, key = lambda k: len(temp4[k]["features"]), reverse = True):
#    sorted_temp4 += [k]
#
## dictionary with compound number as key and associated PC_IDs and features as values
#next_compound_number = max(assigned_compound_numbers + [0])
#to_be_matched = []
#for pcid in sorted_temp4:
#    status = "not matched"
#    if len(temp4[pcid]["compound_number"]) == 1:
#        compound = temp4[pcid]["compound_number"][0]
#        if compound not in temp5:
#            temp5[compound] = {
#                "PC_ID": [pcid],
#                "consensus_features": temp4[pcid]["features"],
#                "features": temp4[pcid]["features"],
#                "average_TIC": temp4[pcid]["average_TIC"]
#            }
#        else:
#            temp5[compound]["PC_ID"] += [pcid]
#            if len(temp5[compound]["consensus_features"]) == 0:
#                temp5[compound]["consensus_features"] = temp4[pcid]["features"]
#            else:
#                temp5[compound]["consensus_features"] = list(set(temp5[compound]["consensus_features"]) & \
#                    set(temp4[pcid]["features"]))
#            for idx, feature in enumerate(temp4[pcid]["features"]):
#                if feature in temp5[compound]["features"]:
#                    temp5[compound]["average_TIC"][temp5[compound]["features"].index(feature)] += \
#                        temp4[pcid]["average_TIC"][idx]
#                else:
#                    temp5[compound]["features"] += [feature]
#                    temp5[compound]["average_TIC"] += [temp4[pcid]["average_TIC"][idx]]
#    else:
#        to_be_matched += [pcid]
#for pcid in to_be_matched:
#    status = "not matched"
#    for compound in temp5.keys():
#        if cosine_score(temp4[pcid], temp5[compound], abundance="average_TIC", features="features") > min_cos_score:
#            temp5[compound]["PC_ID"] += [pcid]
#            if len(temp5[compound]["consensus_features"]) == 0:
#                temp5[compound]["consensus_features"] = temp4[pcid]["features"]
#            else:
#                temp5[compound]["consensus_features"] = list(set(temp5[compound]["consensus_features"]) & \
#                    set(temp4[pcid]["features"]))
#            for idx, feature in enumerate(temp4[pcid]["features"]):
#                if feature in temp5[compound]["features"]:
#                    temp5[compound]["average_TIC"][temp5[compound]["features"].index(feature)] += \
#                        temp4[pcid]["average_TIC"][idx]
#                else:
#                    temp5[compound]["features"] += [feature]
#                    temp5[compound]["average_TIC"] += [temp4[pcid]["average_TIC"][idx]]
#            status = "matched"
#            break
#    if status == "not matched":
#        next_compound_number += 1
#        temp5[next_compound_number] = {
#            "PC_ID": [pcid],
#            "consensus_features": temp4[pcid]["features"],
#            "features": temp4[pcid]["features"],
#            "average_TIC": temp4[pcid]["average_TIC"]
#        }
#
#for compound in temp5:
#    for PC_ID in temp5[compound]["PC_ID"]:
#        compound_table.append({
#            "compound_number": compound,
#            "PC_ID": PC_ID,
#            "features": temp4[PC_ID]["features"] 
#        })



# need another table in database that keeps track of compound numbers so that they stay consistant
#def build_compound_feature(compound_table, pc_id):
#    compound_feature = []
#    for row in pc_id:
#        for compound in compound_table:
#            if row["feature"] in compound["features"] and \
#            row["PC_ID"] == compound["PC_ID"]:
#                for ftrd_id in row["ftrd_id"]:
#                    compound_feature.append({
#                        "ftrd_id": ftrd_id,
#                        "feature_number": row["feature"],
#                        "compound_number": compound["compound_number"]
#                    })
#    
#    return compound_feature


