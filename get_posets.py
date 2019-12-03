import min_interval_posets.posets as mipp
import json,random,sys


def get_samples(names,n,m):
    # randomly pick m sets of n genes from a list, no repeats
    return [tuple(random.sample(names,n)) for _ in range(m)]


def get_single_poset(names,orders):
    # given one set of genes, find the associated partial order
    all_nodes = []
    all_edges = []
    for name in names:
        nodes,edges = orders[name]
        N = len(all_nodes)
        all_nodes.extend(nodes)
        all_edges.extend([(i + N, j + N) for (i, j) in edges])
    extrema,edges = mipp.get_poset(all_nodes, all_edges)
    return extrema,list(edges)


def get_sample_posets(orders,samples):
    # given a list of sets of genes, find all the partial orders and arrange them in a nested dictionary:
    # strain keys epsilon keys set of names keys poset
    posets={}
    for strain,epses in orders.items():
        posets[strain] = {}
        for eps,ords in epses.items():
            posets[strain][eps] = {}
            for names in samples:
                posets[strain][eps][", ".join(names)] = get_single_poset(names,ords)
    return posets


def filter_orders(fname,orders,num_extrema=5):
    # take out all genes that have more than num_extrema in their total orders
    names = [name.split("\n")[0] for name in open(fname)]
    remove = set([])
    for name in names:
        flag = 0
        for strain, epses in orders.items():
            if flag:
                break
            for eps,ords in epses.items():
                if len(ords[name][0]) > num_extrema:
                    remove.add(name)
                    flag = 1
                    break
    return list(set(names).difference(remove))


def filtered_names(inphasefile="3D7_allstrains_outofphase_0.05_two.txt", inorderfile= "inphase_gene_total_orders.json",num_extrema=5,storegenesfile=None):
    # load files and filter to remove genes with too many extrema
    inorders = json.load(open(inorderfile))
    inphase_names = filter_orders(inphasefile,inorders,num_extrema)
    if storegenesfile is not None:
        with open(storegenesfile,"w") as f:
            f.write("\n".join(inphase_names))
    return inorders,inphase_names


def inphase_sample_posets(inorders,inphase_names,num_genes=4,num_samples=1000,outfname="inphase_gene_sample_posets.json"):
    # get the partial orders for num_samples sets of num_genes each
    # save output to a file as well as returning it
    inphase_samples = get_samples(inphase_names,num_genes,num_samples)
    inposets = get_sample_posets(inorders,inphase_samples)
    json.dump(inposets,open(outfname,"w"))
    return inphase_samples, inposets


def permuted_sample_posets(inorders,inphase_names,inphase_samples,inposets,ref="3D7",strains=["D6","FVO","SA250"],outfname="permuted_gene_sample_posets.json",map_outfname="map_names2permutednames.json"):
    # take the information from the inphase partial orders and permute the names to make a base case
    # save the results and the random gene name mapping
    permutenames = list(inphase_names)
    random.shuffle(permutenames)
    permutednamesdict = dict(zip(inphase_names,permutenames))
    backwards = dict(zip(permutenames,inphase_names))
    random_samples = [tuple(permutednamesdict[s] for s in samp) for samp in inphase_samples]
    randomposets = get_sample_posets(inorders,random_samples)
    permutedposets = {}
    permutedposets[ref] = dict(inposets[ref])
    for strain in strains:
        permutedposets[strain] = {}
        for e,d in randomposets[strain].items():
            permutedposets[strain][e] = {}
            for name_string,randposet in d.items():
                names = name_string.split(", ")
                pernames = ", ".join([backwards[name] for name in names])
                new_extrema = []
                for name,extremum in randposet[0]:
                    new_extrema.append([backwards[name],extremum])
                permutedposets[strain][e][pernames] = (tuple(new_extrema),randposet[1])
    json.dump(permutedposets, open(outfname,"w"))
    json.dump(permutednamesdict,open(map_outfname,"w"))


if __name__ == "__main__":
    inorders,inphase_names = filtered_names(inphasefile="3D7_allstrains_outofphase_0.05_two.txt", inorderfile= "inphase_gene_total_orders.json",num_extrema=5)
    inphase_samples,inposets = inphase_sample_posets(inorders,inphase_names,num_genes=4,num_samples=5,outfname="inphase_gene_sample_posets.json")
    permuted_sample_posets(inorders, inphase_names, inphase_samples,inposets,outfname="permuted_gene_sample_posets.json",map_outfname="map_names2permutednames.json")


