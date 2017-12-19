import imp

ct = imp.load_source('ct', './/code//build_compound_table.py')

# if ftrd is correct in new database (it's not currently)
#ftrd = ct.get_ftrd()

# if not correct in new database, load file from old
import pymysql

mydb = pymysql.connect(host='mysql.chpc.utah.edu',
                       user='u6009010',
                       password='3UaUhf7a',
                       db='inga_2015_06_01',
                       charset='utf8mb4',
                       cursorclass=pymysql.cursors.DictCursor)

with mydb.cursor() as cursor:
    sql = "SELECT * FROM `feature_table_raw_data_C18`"
    cursor.execute(sql)
    result = cursor.fetchall()

ftrd = []
for row in result:
    ftrd.append({
        "RT": float(row["RT"]),
        "PC_ID": row["PC_ID"],
        "TIC": float(row["TIC"]),
        "species_code_sample": row["Species_code_sample"],
        "id": row["id"],
        "MZ": float(row["MZ"]),
        "sample": row["sample"]
    })

# build feature table
feature_table = ct.build_feature_table(ftrd, mz_error=0.01, rt_error=0.3)

# write csv of feature table
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//samples_with_polar_features_mzrt.csv", 'w') as file:
    file.write("feature_number,sample,rt,mz,PC_ID,TIC\n")
    for feature in feature_table:
        for sample in feature_table[feature]["sample"]:
            file.write("%d,%s,%f,%f,%s,%f\n" % (feature, sample,\
             feature_table[feature]["rt"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["mz"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["PC_IDs"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["TIC"][feature_table[feature]["sample"].index(sample)]))  

pc_id_table = ct.build_pc_id(feature_table)

compound_table = ct.build_compound_table(pc_id_table, min_cos_score=0.5)

# write csv of compound table
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//polar_compound_feature_table.csv", "w") as file:
    file.write("compound_number,feature_number,TIC,mz,rt\n")
    for compound in compound_table:
        for feature in compound_table[compound]["features"]:
            file.write("%d,%d,%f,%f,%f\n" % (compound, feature, \
            compound_table[compound]["TICs"][compound_table[compound]["features"].index(feature)], \
            feature_table[feature]["avg_mz"], feature_table[feature]["avg_rt"]))

# or read in compound table if it already exists:
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//polar_compound_feature_table.csv", "r") as file:
    compound_table_temp = file.readlines()

compound_table = {}
for row in compound_table_temp:
    compound = row.split(",")[0]
    if compound == 'compound_number':
        continue
    else:
        compound = int(compound)
        if compound in compound_table:
            compound_table[compound]["features"] += [int(row.split(",")[1])]
            compound_table[compound]["TICs"] += [float(row.split(",")[2])]
            compound_table[compound]["avg_rt"] += [float(row.split(",")[4].split("\n")[0].replace('"',''))]
            compound_table[compound]["avg_mz"] += [float(row.split(",")[3])]
        else:
            compound_table[compound] = {
                "features": [int(row.split(",")[1])],
                "TICs": [float(row.split(",")[2])],
                "avg_rt": [float(row.split(",")[4].split("\n")[0].replace('"',''))],
                "avg_mz": [float(row.split(",")[3])]
            }

# load filled features from R Code:
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//samples_with_polar_features_mzrt.csv", "r") as file:
    filled_features_temp = file.readlines()

filled_features = {}
for row in filled_features_temp:
    sample = row.split(",")[4].split("\n")[0].replace('"','')
    if sample == 'sample_name':
        continue
    else:
        if sample in filled_features:
            if int(row.split(",")[0]) not in filled_features[sample]["feature_number"]:
                filled_features[sample]["feature_number"] += [int(row.split(",")[0])]
                filled_features[sample]["TIC"] += [float(row.split(",")[1])]
                filled_features[sample]["actual_rt"] += [float(row.split(",")[3])]
        else:
            filled_features[sample] = {
                "feature_number": [int(row.split(",")[0])],
                "TIC": [float(row.split(",")[1])],
                "actual_rt": [float(row.split(",")[3])]
            }

## calculate percent of TIC for features in each compound
#for compound in compound_table:
#    compound_table[compound]["feature_pcts"] = [x / sum(compound_table[compound]["TICs"]) \
#     for x in compound_table[compound]["TICs"]]
#    compound_table[compound]["rel_feature_abund"] = [x / max(compound_table[compound]["feature_pcts"]) \
#     for x in compound_table[compound]["feature_pcts"]]

filled_comps = ct.fill_compounds(filled_features, compound_table)

# write csv of final compound table
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//polar_compound_table_2017_12_15.csv", "w") as file:
    file.write("compound_sample,compound_number,TIC\n")
    for sample in filled_comps:
        for i, compound in enumerate(filled_comps[sample]["compound"]):
            file.write("%s,%d,%f\n" % (sample, compound, \
            filled_comps[sample]["TIC"][i]))
