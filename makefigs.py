import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
fontsize=20
rc('axes', labelsize=fontsize)    # fontsize of the x and y labels
rc('xtick', labelsize=fontsize)    # fontsize of the tick labels
rc('ytick', labelsize=fontsize)    # fontsize of the tick labels
rc('legend', fontsize=fontsize)    # legend fontsize
from math import sqrt
import json,os, sys, ast
import numpy as np


def get_mean_std(epsilons,normdists):
    msim = []
    ssim = []
    for i in range(len(epsilons)):
        l = [1-n[i] for n in normdists]
        mu = sum(l)/len(l)
        msim.append(mu)
        std = sqrt(sum((a-mu)**2 for a in l) / (len(l)-1))
        ssim.append(std)
    return msim,ssim


def plot_sim_many(epsilons,innormdists,outnormdists,strain,outputfolder,ref):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # ax.set_axisbelow(True)
    in_msim, in_ssim = get_mean_std(epsilons,innormdists)
    out_msim, out_ssim = get_mean_std(epsilons,outnormdists)
    ax.fill_between(epsilons, np.array(in_msim)+np.array(in_ssim), np.array(in_msim)-np.array(in_ssim), facecolor="k", alpha=0.5,interpolate=True)
    plt.plot(epsilons,in_msim,linewidth=3,color="k",zorder=500,label="In phase")
    ax.fill_between(epsilons, np.array(out_msim)+np.array(out_ssim), np.array(out_msim)-np.array(out_ssim), facecolor="r", alpha=0.5,interpolate=True)
    plt.plot(epsilons,out_msim,linewidth=3,color="r",zorder=500,label="Base comparison")
    plt.ylim([0,1.0])
    plt.xlabel(r"$\mathrm{\epsilon}$",fontsize=30)
    plt.ylabel("Similarity")
    plt.legend(loc='lower left')
    ax.title.set_fontsize(30)
    plt.title(strain)
    plt.savefig(os.path.join(outputfolder,"similarity_mean_of_many_{}to{}.pdf".format(strain,ref)), bbox_inches='tight')
    # plt.show()


def plot_results(outputfolder,infname="inphase_gene_sample_pairwise_norm_dists.json",basefname="permuted_gene_sample_pairwise_norm_dists.json",strains=["D6","FVO","SA250"],ref="3D7"):
    inresults = json.load(open(infname))
    epsilons = [float(e) for e in inresults["epsilons"]]
    inresults.pop("epsilons")
    baseresults = json.load(open(basefname))
    baseresults.pop("epsilons")
    in_swapped = {strain : [] for strain in strains}
    for sample,rdict in inresults.items():
        for strain,dists in rdict.items():
            in_swapped[strain].append(dists)
    base_swapped = {strain : [] for strain in strains}
    for sample,rdict in baseresults.items():
        for strain,dists in rdict.items():
            base_swapped[strain].append(dists)
    for strain in strains:
        innormdists = in_swapped[strain]
        outnormdists = base_swapped[strain]
        plot_sim_many(epsilons,innormdists,outnormdists,strain,outputfolder,ref)



if __name__ == "__main__":
    # plot_results("")
    strains = ast.literal_eval(sys.argv[1]) # ["lung","kidney"] for mouse,  ["D6","FVO","SA250"] for malaria
    ref = sys.argv[2] # "liver" for mouse, 3D7 for malaria

    plot_results("",strains=strains,ref=ref)