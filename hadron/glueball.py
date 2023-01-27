
import pandas as pd
import numpy as np
import os
import yaml
import argparse
import glob


parser = argparse.ArgumentParser(
    prog='glueball',
    description="""format glueball correlators and interpolator files in the "hadron" format""")

parser.add_argument('-f', '--inputfile', help="Path to the same yaml input file used for the run")
args = parser.parse_args()

# loading and parsing the input file for the needed info
nd = yaml.load(open(str(args.inputfile)), Loader=yaml.Loader)

T = nd["geometry"]["T"]
nd_omeas = nd["omeas"]
resdir = nd_omeas["res_dir"]
nd_glueball = nd_omeas["glueball"]

save_corr = nd_glueball["correlator"]
nd_interpolator = nd_glueball["interpolator"]

do_APE_smearing = bool(nd_glueball["do_APE_smearing"])
n_APE = []
alpha_APE = 0.0
if (do_APE_smearing):
    nd_APE_smearing = nd_glueball["APE_smearing"]
    n_APE = nd_APE_smearing["n"]
    alpha = nd_APE_smearing["alpha"]
####

rmin = int(nd_interpolator["rmin"])
rmax = int(nd_interpolator["rmax"])
interp_type = nd_interpolator["type"]
interp_save = nd_interpolator["save"]


def hadronize(operator, name_obj, name_ij):
    d1 = resdir+"/"+"glueball/"
    d2 = d1 + operator + "/" + interp_type + "/"
    d2 += "smearAPEn{n_APE}alpha{alpha}/".format(n_APE=n_APE, alpha=alpha)+"/"
    d2 += name_ij+"/"
    print(d2)
    # list all files in d2 and extract configuration numbers from names
    file_names = glob.glob(d2+"/"+name_obj+"_*")
    file_index = [int(x.split(name_obj+"_")[1]) for x in file_names]
    #
    # load .dat file with indices of formatted configurations, if existing
    iconf_file = d2 + "/iconfs.dat"
    # already formatted configurations
    file_index_old = []
    if os.path.exists(iconf_file):
        print("Old configurations found")
        file_index_old = pd.read_csv(
            iconf_file, sep=" ", header=None)[0].tolist()
    ####
    n_conf_old = len(file_index_old)
    # new matrices for a each J^{PC}
    # dict({"pp": pd.DataFrame(), "pm": pd.DataFrame(), "mp": pd.DataFrame(), "mm":pd.DataFrame()})
    df_new = dict()
    print("Reading old data if present")
    for j_pc in ["pp", "pm", "mp", "mm"]:
        # loop over all J^{PC} quantum numbers
        df_new[j_pc] = pd.DataFrame(np.zeros(shape=(1, T+1)), columns=["i"]+[
            str(t) for t in range(T)], index=["remove"])
        df_new[j_pc] = df_new[j_pc].astype({"i": int})
        if n_conf_old > 0:
            path_old = d2+"/"+j_pc+".dat"
            if os.path.exists(iconf_file) and (not os.path.exists(path_old)):
                print("Error! Something bad happened, please redo the online measurements:\n",
                      iconf_file, "\nexists, but\n", path_old, "\ndoesn't.")
                raise ValueError
            ####
            df_old = pd.read_csv(path_old, sep=" ", dtype={"i": int})
            df_new[j_pc] = pd.concat([df_new[j_pc], df_old], axis=0)
            # print(df_new[j_pc])
        ####
    ####
    print("loop over new and updated files")
    for i_f in range(len(file_index)):
        path_i = file_names[i_f]
        print(i_f, path_i)
        idx_i = file_index[i_f]
        df_i = pd.read_csv(path_i, sep=" ")
        for j_pc in ["pp", "pm", "mp", "mm"]:
            # loop over all J^{PC} quantum numbers
            arr1 = df_i[j_pc].to_numpy()
            if idx_i in file_index_old:
                # removing old configuration
                df_new[j_pc] = df_new[j_pc][df_new[j_pc]["i"] != idx_i]
            ####
            df1 = pd.DataFrame(arr1).transpose()
            df1.insert(0, "i", idx_i)
            df1.columns = ["i"] + [str(t) for t in range(T)]
            # update new array with the new data
            df_new[j_pc] = pd.concat([df_new[j_pc], df1])
            ####
        ####
    ####
    print("Cleaning up the dataframe format before exporting")
    save_iconf = True 
    for j_pc in ["pp", "pm", "mp", "mm"]:
        df_new[j_pc] = df_new[j_pc].drop("remove")
        save_df_new = (df_new[j_pc].shape[0]>0)
        save_iconf = save_iconf and save_df_new
        # saving
        if save_df_new:
            df_new[j_pc] = df_new[j_pc].sort_values(by=['i'])
            df_new[j_pc].to_csv(d2 + "/"+j_pc+".dat", sep=" ", index=False)
    ####
    print("Saving the new iconf file")
    if save_iconf:
        pd.DataFrame(sorted(file_index)).to_csv(
            d2+"iconfs.dat", header=False, sep=" ", index=False)
    ####
    print("Removing data in the old format")
    for p in file_names:
        print("removing:", p)
        os.remove(p)
    ####
####


for i in range(rmin, rmax+1):
    hadronize("interpolator", "phi", str(i))
    for j in range(rmin, i+1):
        hadronize("correlator", "C_glueball", str(i)+"_"+str(j))

