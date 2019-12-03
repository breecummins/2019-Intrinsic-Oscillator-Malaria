import os,json


def common_params(outputfolder):
    param_dict = {}
    param_dict["gto_outfname"] = os.path.join(outputfolder,"inphase_gene_total_orders.json")
    param_dict["genes_used"] = os.path.join(outputfolder,"inphase_genes_used.txt")
    param_dict["epsilons"] = [n/100 for n in range(6,11,1)]
    param_dict["num_extrema"] = 8
    param_dict["num_genes"]=6
    param_dict["num_samples"]=5000
    param_dict["gp_inphase_outfname"]=os.path.join(outputfolder,"inphase_gene_sample_posets.json")
    param_dict["gp_base_outfname"]=os.path.join(outputfolder,"permuted_gene_sample_posets.json")
    param_dict["gp_map_outfname"]=os.path.join(outputfolder,"map_names2permutednames.json")
    param_dict["gpd_inphase_outfname"]=os.path.join(outputfolder,"inphase_gene_sample_pairwise_norm_dists.json")
    param_dict["gpd_permute_outfname"]=os.path.join(outputfolder,"permuted_gene_sample_pairwise_norm_dists.json")
    return param_dict


def define_params_mal(outputfolder):
    # these parameters will be saved with the output of each run, so this function can be overwritten as needed
    param_dict = common_params(outputfolder)
    param_dict["inphasefile"]="3D7_allstrains_outofphase_0.05_two.txt"
    param_dict["strainfiles"] = {"3D7": "downsampled_3D7_expression_wrapped_offset.csv","D6" : "downsampled_D6_expression_wrapped_offset.csv", "FVO" : "downsampled_FVO_expression_wrapped_offset.csv", "SA250":"downsampled_SA250_expression_wrapped_offset.csv"}
    param_dict["reference"] = "3D7"
    pfname = os.path.join(outputfolder,"params.json")
    json.dump(param_dict,open(pfname,"w"))
    return outputfolder,pfname


def define_params_mouse(outputfolder):
    # these parameters will be saved with the output of each run, so this function can be overwritten as needed
    param_dict = common_params(outputfolder)
    param_dict["inphasefile"] = "mouse_gene_in_phase.txt"
    param_dict["strainfiles"] = {"lung" : "lung_in_phase_0_05.tsv","liver" : "liver_in_phase_0_05.tsv","kidney" : "kidney_in_phase_0_05.tsv"}
    param_dict["reference"] = "liver"
    pfname = os.path.join(outputfolder,"params.json")
    json.dump(param_dict,open(pfname,"w"))
    return outputfolder,pfname


