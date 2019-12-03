import min_interval_posets.curve as mipc
import min_interval_posets.posets as mipp
import pandas,json,sys


def extractdata(filename):
    # Turn time series file into pandas dataframe
    file_type = filename.split(".")[-1]
    if file_type == "tsv":
        df = pandas.read_csv(open(filename),delim_whitespace=True, dtype=str)
    elif file_type == "csv":
        df = pandas.read_csv(open(filename), dtype=str)
    else:
        raise ValueError("File type not recognized. Require .tsv or .csv.")
    return list(df)[1:],df.values


def row(filename):
    # parse row formatted time series file (one gene per row)
    # create a dictionary keying gene name to Curve object from min_interval_posets
    times,data = extractdata(filename)
    times = [float(n) for n in times]
    names = data[:,0]
    curves = [mipc.Curve(data[k,1:].astype(float),times,True) for k in range(data.shape[0])]
    return dict(zip(names,curves))


def col(filename):
    # parse column formatted time series file (one gene per column)
    # create a dictionary keying gene name to Curve object from min_interval_posets
    names,data = extractdata(filename)
    times = data[:,0]
    curves = [mipc.Curve(data[:,k].astype(float),times,True) for k in range(1,data.shape[1])]
    return dict(zip(names,curves))


def make_curves(inphasefile="3D7_allstrains_outofphase_0.05_two.txt",strainfiles = {"3D7":
            "3D7_expression_wrapped_offset.csv","D6" : "D6_expression_wrapped_offset.csv", "FVO" :
            "FVO_expression_wrapped_offset.csv", "SA250":"SA250_expression_wrapped_offset.csv"},func=row):
    # Form a nested dictionary first keyed by strain and then by gene name to a Curve object from min_interval_posets
    # The gene names are limited to those in the inphasefile
    # The keyword "func" specifies row or column format of the time series files in strainfiles.
    # All files must all have the same format.
    inphase = [name.split("\n")[0] for name in open(inphasefile)]
    incurves = {}
    for strain,sf in strainfiles.items():
        all_curves = func(sf)
        incurves[strain] = {}
        for name in inphase:
            incurves[strain][name] = all_curves[name]
    return incurves


def make_orders(curves,epsilons):
    # takes the nested dictionary output of make_curves() and a list of epsilons between 0 and 0.5
    # returns a nested dictionary of strain keying epsilon keying name keying Curve object from min_interval_posets
    orders = {}
    for strain in curves.keys():
        orders[strain] = {}
        for eps in epsilons:
            orders[strain][eps] = {}
            for name,curve in curves[strain].items():
                orders[strain][eps][name] = mipp.get_total_order(name, curve, eps)
    return orders


def get_total_orders(epsilons,incurves,outfname="inphase_gene_total_orders.json"):
    # takes the output of make_curves()
    # saves a nested dictionary of strain keying epsilon keying name keying Curve object from min_interval_posets
    json.dump(make_orders(incurves,epsilons),open(outfname,"w"))


if __name__ == "__main__":
    incurves = make_curves()
    epsilons = [n/100 for n in range(5,10,1)]
    get_total_orders(epsilons,incurves,outfname="inphase_gene_total_orders.json")
