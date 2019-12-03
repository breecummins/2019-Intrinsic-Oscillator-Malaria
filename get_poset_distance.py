import min_interval_posets.poset_distance as pdist
import json,sys
from multiprocessing import Pool
from functools import partial


def pos_dist(posets,N=4,ref="3D7"):
    # set up parallel computation to find pairwise distances between reference strain and other strains
    strains = list(posets.keys())
    strains.remove(ref)
    epsilons = sorted(list(posets[strains[0]].keys()))
    samples = list(posets[strains[0]][epsilons[0]].keys())
    pool = Pool(N)
    output = pool.map(partial(inner,epsilons,strains,posets,ref),samples)
    results = dict(zip(samples,output))
    results.update({"epsilons" : epsilons})
    return results


def inner(epsilons,strains,posets,ref,samp):
    # helper function for the parallel process
    results = dict(zip(strains,[[] for _ in range(len(strains))]))
    for eps in epsilons:
        refgraph = pdist.poset_to_nx_graph(posets[ref][eps][samp])
        for strain in strains:
            print(samp, eps, strain) #keep track of computation progress
            straingraph = pdist.poset_to_nx_graph(posets[strain][eps][samp])
            results[strain].append(pdist.normalized_dag_distance(refgraph, straingraph))
    return results


# def in_and_out_dists(inposfile="inphase_gene_sample_posets.json",outposfile="outphase_gene_sample_posets.json"):
#     # OBSOLETE
#     inposets = json.load(open(inposfile))
#     outposets = json.load(open(outposfile))
#     print("In phase\n")
#     indists = pos_dist(inposets)
#     print("Out of phase\n")
#     outdists = pos_dist(outposets)
#     json.dump(indists,open("inphase_gene_sample_pairwise_norm_dists.json","w"))
#     json.dump(outdists, open("outphase_gene_sample_pairwise_norm_dists.json","w"))
#
#
# def in_and_mixed_dists(inposfile="inphase_gene_sample_posets.json",mixedposfile="mixedphase_gene_sample_posets.json"):
#     # OBSOLETE
#     inposets = json.load(open(inposfile))
#     mixedposets = json.load(open(mixedposfile))
#     print("In phase\n")
#     indists = pos_dist(inposets)
#     print("Mixed phase\n")
#     mixeddists = pos_dist(mixedposets)
#     json.dump(indists,open("inphase_gene_sample_pairwise_norm_dists.json","w"))
#     json.dump(mixeddists, open("mixedphase_gene_sample_pairwise_norm_dists.json","w"))

def in_and_permuted_dists(N=4,ref="3D7",inposfile="inphase_gene_sample_posets.json",permutedfile="permuted_gene_sample_posets.json",inphase_outfname="inphase_gene_sample_pairwise_norm_dists.json",permute_outfname="permuted_gene_sample_pairwise_norm_dists.json"):
    # find distances for in phase + base case
    # names are currently set for the permuted case, but one could make the output file name variable for different base cases
    inposets = json.load(open(inposfile))
    indists = pos_dist(inposets,N,ref)
    perposets = json.load(open(permutedfile))
    perdists = pos_dist(perposets,N,ref)
    json.dump(indists,open(inphase_outfname,"w"))
    json.dump(perdists,open(permute_outfname,"w"))


if __name__ == "__main__":
    if len(sys.argv) > 1:
        in_and_permuted_dists(sys.argv[1],inposfile="inphase_gene_sample_posets.json",
                              permutedfile="permuted_gene_sample_posets.json",
                              inphase_outfname="inphase_gene_sample_pairwise_norm_dists.json",
                              permute_outfname="permuted_gene_sample_pairwise_norm_dists.json")
    else:
        in_and_permuted_dists(inposfile="inphase_gene_sample_posets.json",permutedfile="permuted_gene_sample_posets.json",inphase_outfname="inphase_gene_sample_pairwise_norm_dists.json",permute_outfname="permuted_gene_sample_pairwise_norm_dists.json")
