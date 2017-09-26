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
def match_features(raw_data, mz_error, rt_error):
   feature_table = {}
   current_feature_number = 0
   for datum in raw_data:
        pcid = datum["PC_ID"]
        mz = datum["MZ"]
        rt = datum["RT"]
        status = "not matched"
        for feature in feature_table.keys():
            if abs(feature_table[feature]["avg_mz"] - mz) <= mz_error and \
                abs(feature_table[feature]["avg_rt"] - rt) <= rt_error:
                feature_table[feature]["mz"] += [mz]
                feature_table[feature]["rt"] += [rt]
                feature_table[feature]["avg_mz"] = sum(feature_table[feature]["mz"])/len(feature_table[feature]["mz"])
                feature_table[feature]["avg_rt"] = sum(feature_table[feature]["rt"])/len(feature_table[feature]["rt"])
                status = "matched"
                break
        if status == "not matched":
            current_feature_number += 1
            feature_table[current_feature_number] = {
                "mz": [mz],
                "rt": [rt],
                "avg_mz": mz,
                "avg_rt": rt,
                }


# build pc_id table 
def build_pc_id(raw_data, mz_error, rt_error, compound_feature = "missing"):
    temp = {}
    temp2 = {}
    feature_table = {}
    temp3 = {}
    pc_id_1 = []
    pc_id_2 = []

    if compound_feature == "missing":
        compound_feature = [{'feature_number': None, 'ftrd_id': None, 'compound_number': None}]

    temp7 = {}
    for row in raw_data:
        temp7[row["id"]] = row
    temp8 = {}
    for row in compound_feature:
        temp8[row["ftrd_id"]] = row
    raw_data_2 = []
    for datum in temp7:
        if datum in temp8:
            raw_data_2.append({
                "RT": temp7[datum]["RT"],
                "MZ": temp7[datum]["MZ"],
                "sample": temp7[datum]["sample"],
                "PC_ID": temp7[datum]["PC_ID"],
                "TIC": temp7[datum]["TIC"],
                "species_code_sample": temp7[datum]["species_code_sample"],
                "id": datum,
                "feature_number": [temp8[datum]["feature_number"]],
                "compound_number": [temp8[datum]["compound_number"]]
            })
        else:
            raw_data_2.append({
                "RT": temp7[datum]["RT"],
                "MZ": temp7[datum]["MZ"],
                "sample": temp7[datum]["sample"],
                "PC_ID": temp7[datum]["PC_ID"],
                "TIC": temp7[datum]["TIC"],
                "species_code_sample": temp7[datum]["species_code_sample"],
                "id": datum,
                "feature_number": [None],
                "compound_number": [None]
            })

    for datum in raw_data_2:
        pcid_mzrt = str(datum["PC_ID"]) + '_' + str(datum["MZ"]) + '_' + str(datum["RT"])
        if pcid_mzrt in temp:
            temp[pcid_mzrt]["TIC"] += [datum["TIC"]]
            temp[pcid_mzrt]["row_ids"] += [datum["id"]]
            temp[pcid_mzrt]["feature_number"] = list(set(temp[pcid_mzrt]["feature_number"]) | \
             set(datum["feature_number"]))
            temp[pcid_mzrt]["compound_number"] = list(set(temp[pcid_mzrt]["compound_number"]) | \
             set(datum["compound_number"]))
        else:
            temp[pcid_mzrt] = {
                "TIC": [datum["TIC"]],
                "row_ids": [datum["id"]],
                "feature_number": datum["feature_number"],
                "compound_number": datum["compound_number"]
            }
        if datum["PC_ID"] in temp2:
            temp2[datum["PC_ID"]]["Total_TIC"] += datum["TIC"]
        else:
            temp2[datum["PC_ID"]] = {
                "Total_TIC": datum["TIC"]
            }

    current_feature_number = max([x["feature_number"] for x in compound_feature] + [0])
    to_be_matched = []
    for pcidmzrt in temp.keys():
        pcid = '_'.join([pcidmzrt.split('_')[x] for x in [0,1]])
        mz = float(pcidmzrt.split('_')[2])
        rt = float(pcidmzrt.split('_')[3])
        status = "not matched"
        feature = temp[pcidmzrt]["feature_number"][0]
        if  feature is not None:
            if feature in feature_table.keys():
                if pcid in feature_table[feature]["PC_ID"]:
                    feature_table[feature][pcid]["ftrd_ids"] += temp[pcidmzrt]["row_ids"]
                    feature_table[feature][pcid]["TIC"] += temp[pcidmzrt]["TIC"]
                    feature_table[feature][pcid]["compound_number"] = \
                     list(set(feature_table[feature][pcid]["compound_number"]) | \
                     set(temp[pcidmzrt]["compound_number"]))
                else:
                    feature_table[feature]["mz"] += [mz]
                    feature_table[feature]["rt"] += [rt]
                    feature_table[feature]["PC_ID"] += [pcid]
                    feature_table[feature][pcid] = {
                        "ftrd_ids": temp[pcidmzrt]["row_ids"],
                        "TIC": temp[pcidmzrt]["TIC"],
                        "compound_number": temp[pcidmzrt]["compound_number"]
                    }
            else:
                feature_table[feature] = {
                    "mz": [mz],
                    "rt": [rt],
                    "PC_ID": [pcid],
                    pcid: {
                        "ftrd_ids": temp[pcidmzrt]["row_ids"],
                        "TIC": temp[pcidmzrt]["TIC"],
                        "compound_number": temp[pcidmzrt]["compound_number"]
                    }
                }
        else:
            to_be_matched += [pcidmzrt]
    for pcidmzrt in to_be_matched:
        pcid = '_'.join([pcidmzrt.split('_')[x] for x in [0,1]])
        mz = float(pcidmzrt.split('_')[2])
        rt = float(pcidmzrt.split('_')[3])
        status = "not matched"
        for feature in feature_table.keys():
            if abs(sum(feature_table[feature]["mz"])/len(feature_table[feature]["mz"]) - mz) <= mz_error and \
                abs(sum(feature_table[feature]["rt"])/len(feature_table[feature]["rt"]) - rt) <= rt_error:
                if pcid in feature_table[feature]["PC_ID"]:
                    feature_table[feature][pcid]["ftrd_ids"] += temp[pcidmzrt]["row_ids"]
                    feature_table[feature][pcid]["TIC"] += temp[pcidmzrt]["TIC"]
                else:
                    feature_table[feature]["mz"] += [mz]
                    feature_table[feature]["rt"] += [rt]
                    feature_table[feature]["PC_ID"] += [pcid]
                    feature_table[feature][pcid] = {
                        "ftrd_ids": temp[pcidmzrt]["row_ids"],
                        "TIC": temp[pcidmzrt]["TIC"],
                        "compound_number": []
                    }
                status = "matched"
                break
        if status == "not matched":
            current_feature_number += 1
            feature_table[current_feature_number] = {
                "mz": [mz],
                "rt": [rt],
                "PC_ID": [pcid],
                pcid: {
                    "ftrd_ids": temp[pcidmzrt]["row_ids"],
                    "TIC": temp[pcidmzrt]["TIC"],
                    "compound_number": []
                }
            }

    for feature in feature_table:
        for idx, pcid in enumerate(feature_table[feature]["PC_ID"]):
            pc_id_2.append({
                "PC_ID": pcid,
                "feature": feature,
                "rt": feature_table[feature]["rt"][idx],
                "mz": feature_table[feature]["mz"][idx],
                "average_TIC": sum(feature_table[feature][pcid]["TIC"]) / \
                 len(feature_table[feature][pcid]["TIC"]),
                "percent_TIC": sum(feature_table[feature][pcid]["TIC"]) / temp2[pcid]["Total_TIC"],
                "ftrd_id": feature_table[feature][pcid]["ftrd_ids"],
                "compound_number": feature_table[feature][pcid]["compound_number"]
            })
    
    return pc_id_2

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

def build_compound_table(pc_id, min_cos_score):
    temp4 = {}
    sorted_temp4 = []
    temp5 = {}
    compound_table = []

    assigned_compound_numbers = [x["compound_number"][0] for x in pc_id if len(x["compound_number"]) > 0]

    for feature in pc_id:
        current_pcid = feature["PC_ID"]
        if current_pcid in temp4:
            if feature["feature"] in temp4[current_pcid]["features"]:
                temp4[current_pcid]["average_TIC"][temp4[current_pcid]["features"].index(feature["feature"])] += \
                 feature["average_TIC"]
            else:
                temp4[current_pcid]["features"] += [feature["feature"]]
                temp4[current_pcid]["average_TIC"] += [feature["average_TIC"]]
        else:
            temp4[current_pcid] = {
                "features": [feature["feature"]],
                "average_TIC": [feature["average_TIC"]],
                "compound_number": feature["compound_number"]
            }

    for k in sorted(temp4, key = lambda k: len(temp4[k]["features"]), reverse = True):
        sorted_temp4 += [k]
    
    # dictionary with compound number as key and associated PC_IDs and features as values
    next_compound_number = max(assigned_compound_numbers + [0])
    to_be_matched = []
    for pcid in sorted_temp4:
        status = "not matched"
        if len(temp4[pcid]["compound_number"]) == 1:
            compound = temp4[pcid]["compound_number"][0]
            if compound not in temp5:
                temp5[compound] = {
                    "PC_ID": [pcid],
                    "consensus_features": temp4[pcid]["features"],
                    "features": temp4[pcid]["features"],
                    "average_TIC": temp4[pcid]["average_TIC"]
                }
            else:
                temp5[compound]["PC_ID"] += [pcid]
                if len(temp5[compound]["consensus_features"]) == 0:
                    temp5[compound]["consensus_features"] = temp4[pcid]["features"]
                else:
                    temp5[compound]["consensus_features"] = list(set(temp5[compound]["consensus_features"]) & \
                        set(temp4[pcid]["features"]))
                for idx, feature in enumerate(temp4[pcid]["features"]):
                    if feature in temp5[compound]["features"]:
                        temp5[compound]["average_TIC"][temp5[compound]["features"].index(feature)] += \
                            temp4[pcid]["average_TIC"][idx]
                    else:
                        temp5[compound]["features"] += [feature]
                        temp5[compound]["average_TIC"] += [temp4[pcid]["average_TIC"][idx]]
        else:
            to_be_matched += [pcid]
    for pcid in to_be_matched:
        status = "not matched"
        for compound in temp5.keys():
            if cosine_score(temp4[pcid], temp5[compound], abundance="average_TIC", features="features") > min_cos_score:
                temp5[compound]["PC_ID"] += [pcid]
                if len(temp5[compound]["consensus_features"]) == 0:
                    temp5[compound]["consensus_features"] = temp4[pcid]["features"]
                else:
                    temp5[compound]["consensus_features"] = list(set(temp5[compound]["consensus_features"]) & \
                        set(temp4[pcid]["features"]))
                for idx, feature in enumerate(temp4[pcid]["features"]):
                    if feature in temp5[compound]["features"]:
                        temp5[compound]["average_TIC"][temp5[compound]["features"].index(feature)] += \
                            temp4[pcid]["average_TIC"][idx]
                    else:
                        temp5[compound]["features"] += [feature]
                        temp5[compound]["average_TIC"] += [temp4[pcid]["average_TIC"][idx]]
                status = "matched"
                break
        if status == "not matched":
            next_compound_number += 1
            temp5[next_compound_number] = {
                "PC_ID": [pcid],
                "consensus_features": temp4[pcid]["features"],
                "features": temp4[pcid]["features"],
                "average_TIC": temp4[pcid]["average_TIC"]
            }
    
    for compound in temp5:
        for PC_ID in temp5[compound]["PC_ID"]:
            compound_table.append({
                "compound_number": compound,
                "PC_ID": PC_ID,
                "features": temp4[PC_ID]["features"] 
            })
    
    return compound_table

# need another table in database that keeps track of compound numbers so that they stay consistant
def build_compound_feature(compound_table, pc_id):
    compound_feature = []
    for row in pc_id:
        for compound in compound_table:
            if row["feature"] in compound["features"] and \
            row["PC_ID"] == compound["PC_ID"]:
                for ftrd_id in row["ftrd_id"]:
                    compound_feature.append({
                        "ftrd_id": ftrd_id,
                        "feature_number": row["feature"],
                        "compound_number": compound["compound_number"]
                    })
    
    return compound_feature


