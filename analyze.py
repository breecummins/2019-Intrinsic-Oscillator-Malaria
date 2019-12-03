import json
from math import sqrt

def load_files(infname="inphase_gene_sample_pairwise_norm_dists.json",basefname="permuted_gene_sample_pairwise_norm_dists.json",strains=["D6","FVO","SA250"]):
    inresults = json.load(open(infname))
    epsilons = [float(e) for e in inresults["epsilons"]]
    inresults.pop("epsilons")
    baseresults = json.load(open(basefname))
    baseresults.pop("epsilons")
    inphase = {strain : [] for strain in strains}
    baseline = {strain : [] for strain in strains}
    diff = {strain : [] for strain in strains}
    alpha = {strain : 0 for strain in strains}
    for sample,rdict in inresults.items():
        for strain,dists in rdict.items():
            inphase_sim = [1- d for d in dists]
            baseline_sim = [1-b for b in baseresults[sample][strain]]
            inphase[strain].extend(inphase_sim)
            baseline[strain].extend(baseline_sim)
            diff[strain].extend([(a_i) - (b_i) for a_i, b_i in zip(inphase_sim,baseline_sim) ])
            alpha[strain]+= sum([(a_i) > (b_i) for a_i, b_i in zip(inphase_sim,baseline_sim) ])
    analysis = {strain : {} for strain in strains}
    for strain in strains:
        in_total_mean, in_total_std = get_mean_std(inphase[strain])
        base_total_mean, base_total_std = get_mean_std(baseline[strain])
        diff_mean,diff_std = get_mean_std(diff[strain])
        analysis[strain] = {"total baseline mean" : base_total_mean, "total baseline std": base_total_std,"diff mean" : diff_mean, "diff std" : diff_std, "total inphase mean" : in_total_mean, "total inphase std" : in_total_std, "alpha" : alpha[strain]/len(inphase[strain])}
    for a,r in analysis.items():
        print(a)
        print(r)
        print("\n")


def get_mean_std(l):
    mu = sum(l)/len(l)
    return mu, sqrt(sum((a-mu)**2 for a in l) / (len(l)-1))


if __name__ == "__main__":
    load_files(strains = ["D6","FVO","SA250"])
    # load_files(strains = ["lung","kidney"])