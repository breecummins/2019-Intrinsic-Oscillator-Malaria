import get_total_orders as gto
import get_posets as gp
import get_poset_distance as gpd
import makefigs as mf
import json,subprocess,sys,os
# keep the following here -- they need to appear in globals()
from define_params import define_params_mal, define_params_mouse



def create_output_folder(prefix="results"):
    # datetime stamp the output to avoid namespace collision
    datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'],shell=True).decode(sys.stdout.encoding).strip()
    outputfolder = prefix+datetime
    os.makedirs(outputfolder)
    return outputfolder


def shift_name(filename):
    if type(filename) == dict:
        for strain in filename :
            i = filename[strain].index('.')
            filename[strain] = filename[strain][:i] + '_shifted.csv'
    else:
        i = filename.index('.')
        filename = filename[:i] + '_shifted.csv'
    return filename


def main(prefix,define_params,N=4,ts_format = gto.row, shift = True):
    # runs through get_total_orders.py, get_posets.py, get_poset_distance.py, and makefigs.py in order
    # see modules for comments
    # the basic form is that each step saves a file that feeds into the next step
    # all files and parameters are saved to an output folder in the current directory
    # prefix = results folder prefix
    # N = number of cores for parallel process

    outputfolder,pfname = globals()[define_params](create_output_folder(prefix))
    pd = json.load(open(pfname))
    curves = gto.make_curves(inphasefile=pd["inphasefile"],strainfiles = pd["strainfiles"],func=ts_format)
    shifted_curves = gto.make_curves(inphasefile=pd["inphasefile"],strainfiles = shift_name(pd["strainfiles"]),func=ts_format)
    gto.get_total_orders(pd["epsilons"],curves,outfname=pd["gto_outfname"])
    gto.get_total_orders(pd["epsilons"],shifted_curves,outfname=shift_name(pd["gto_outfname"]))
    print("Total orders complete.")

    inorders,inphase_names = gp.filtered_names(inphasefile=pd["inphasefile"], inorderfile= pd["gto_outfname"],num_extrema=pd["num_extrema"],storegenesfile=pd["genes_used"])
    inorders_shifted,_     = gp.filtered_names(inphasefile=pd["inphasefile"], inorderfile= shift_name(pd["gto_outfname"]),num_extrema=pd["num_extrema"])
    # sys.exit()
    inphase_samples,inposets = gp.inphase_sample_posets(inorders,inphase_names,num_genes=pd["num_genes"],num_samples=pd["num_samples"],outfname=pd["gp_inphase_outfname"])
    # inphase_samples_shifted,inposets_shifted = gp.inphase_sample_posets(inorders_shifted,inphase_names,num_genes=pd["num_genes"],num_samples=pd["num_samples"],outfname=shift_name(pd["gp_inphase_outfname"]))
    strains = list(pd["strainfiles"].keys())
    strains.remove(pd["reference"])
    if shift == True:
        gp.permuted_sample_posets(inorders_shifted, inphase_names, inphase_samples,inposets,ref=pd["reference"],strains=strains,outfname=pd["gp_base_outfname"],map_outfname=pd["gp_map_outfname"])
    else:
        gp.permuted_sample_posets(inorders, inphase_names, inphase_samples,inposets,ref=pd["reference"],strains=strains,outfname=pd["gp_base_outfname"],map_outfname=pd["gp_map_outfname"])
    print("Posets complete.")
    gpd.in_and_permuted_dists(N,ref=pd["reference"],inposfile=pd["gp_inphase_outfname"],permutedfile=pd["gp_base_outfname"],inphase_outfname=pd["gpd_inphase_outfname"],permute_outfname=pd["gpd_permute_outfname"])
    print("Distances complete.")
    mf.plot_results(outputfolder,infname=pd["gpd_inphase_outfname"],
                 basefname=pd["gpd_permute_outfname"], strains=strains)


if __name__ == "__main__":
    # command line call
    # python malaria_scripts.py <define params function> <number of cores>
    main("results",sys.argv[1],int(sys.argv[2]))
