import imp

ct = imp.load_source('ct', './/code//build_compound_table.py')

# if ftrd is correct in new database
ftrd = ct.get_ftrd()

# if not correct in new database, load file from old as in "xfer_ftrd_old_to_new_db.py"
import pymysql
import compound_table as ct

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

ftrd1 = []
for row in result:
    ftrd1.append({
        "RT": float(row["RT"]),
        "PC_ID": row["PC_ID"],
        "TIC": float(row["TIC"]),
        "species_code_sample": row["Species_code_sample"],
        "id": row["id"],
        "MZ": float(row["MZ"]),
        "sample": row["sample"]
    })

# build feature table
feature_table = ct.build_feature_table(ftrd1, mz_error = 0.01, rt_error = 0.3)

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
pc_id_table_orig = pc_id_table

compound_table = ct.build_compound_table(pc_id_table, min_cos_score = 0.5)

with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//polar_compound_feature_table.csv", "w") as file:
    file.write("compound_number,feature_number,TIC,mz,rt\n")
    for compound in compound_table:
        for feature in compound_table[compound]["features"]:
            file.write("%d,%d,%f,%f,%f\n" % (compound, feature, \
            compound_table[compound]["TICs"][compound_table[compound]["features"].index(feature)], \
            feature_table[feature]["avg_mz"], feature_table[feature]["avg_rt"]))

# or read in compound table if it already exists:
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//compound_feature_table.csv", "r") as file:
    compound_table_temp = file.readlines()

compound_table = {}
for row in compound_table_temp:
    compound = row.split(",")[0]
    if compound == 'compound_number':
        next
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
with open("K://GY_LAB_FILES//github_repositories//chem_similarity//results//LA4_filled_features_ppm_2017_11_17.csv", "r") as file:
    filled_features_temp = file.readlines()

filled_features = {}
for row in filled_features_temp:
    sample = row.split(",")[4].split("\n")[0].replace('"','')
    if sample == 'sample_name':
        next
    else:
        if sample in filled_features:
            filled_features[sample]["feature_number"] += [int(row.split(",")[0])]
            filled_features[sample]["TIC"] += [float(row.split(",")[1])]
            filled_features[sample]["actual_rt"] += [float(row.split(",")[3])]
        else:
            filled_features[sample] = {
                "feature_number": [int(row.split(",")[0])],
                "TIC": [float(row.split(",")[1])],
                "actual_rt": [float(row.split(",")[3])]
            }

# calculate percent of TIC for features in each compound
for compound in compound_table:
    compound_table[compound]["feature_pcts"] = [x / sum(compound_table[compound]["TICs"]) \
     for x in compound_table[compound]["TICs"]]
    compound_table[compound]["rel_feature_abund"] = [x / max(compound_table[compound]["feature_pcts"]) \
     for x in compound_table[compound]["feature_pcts"]]

filled_compounds_05 = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.5)
filled_compounds_03 = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)
# changed build_compound_table_code to require sample to contain at least 20% of peaks to have a compound.
filled_comps_03_min02 = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)
# changed parameters for fill features (using 25ppm error for mz, at least 1000 TIC rather than 1500)
filled_comps_ppm = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)
# change parameters for fill compounds -- must contain at least 10% of features to have a compound.
# 20% was weeding out too many that had a compound at low abundance (see compound 502)
filled_comps_ppm_min01 = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)
# change parameters to not remove features after they are used once.
filled_comps = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)
filled_comp_2 = ct.fill_compounds(filled_features, compound_table, min_cos_score = 0.3)

with open("K://GY_LAB_FILES//github_repositories//chem_similarity//data//LA4_filled_compound_table_2017_11_30.csv", "w") as file:
    file.write("compound_sample,compound_number,TIC\n")
    for sample in filled_comp_2:
        for i, compound in enumerate(filled_comp_2[sample]["compound"]):
            file.write("%s,%d,%f\n" % (sample, compound, \
            filled_comp_2[sample]["TIC"][i]))


filled_compounds_1[filled_compounds_1.keys()[0]]

sample = filled_features.keys()[2]
compound = 500
min_cos_score = 0.3
filled_compounds = {}

compound_table[841]
[x for x in compound_table.keys() if 11185 in compound_table[x]["features"]]

[x for x in filled_compounds_05.keys() if 500 in filled_compounds_05[x]["compound"]]
[x for x in filled_compounds_03.keys() if 502 in filled_compounds_03[x]["compound"]]
[x for x in filled_comps_03_min02.keys() if 502 in filled_comps_03_min02[x]["compound"]]
[x for x in filled_comps_ppm.keys() if 1 in filled_comps_ppm[x]["compound"]]
[x for x in filled_comps_ppm_min01.keys() if 898 in filled_comps_ppm_min01[x]["compound"]]
filled_comps_03_min02[filled_comps_03_min02.keys()[0]]
filled_comps_ppm_min01["T76_2163_undiluted"]
compound_table[502]
# compound 1 orginally came from pc group LA17_3
# which is NOT the same compound as is in most of the samples that are showing up as containing compound 1. How was it hijacked?
pc_id_table["LA17_3"]



filled_compounds = filled_compounds_1
filled_compounds['IngU_1933']["compound"]
filled_compounds['IngU_1402_2']["compound"]
list(set(filled_comps_ppm_min01['IngU_1933']["compound"]).intersection(filled_comps_ppm_min01['IngU_1402_2']["compound"]))
compound_table[898]