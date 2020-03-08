import pandas as pd
import seaborn as sns
from hashlib import md5
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import os
import pickle
import numpy as np
from scipy.stats import percentileofscore

# read domains from .ann file and store as pd dataframe
def read_domains(file):
    domains = pd.read_csv(file,sep="\t",header=None)
    domains.rename(columns={0:"chr1",1:"start1",2:"end1",3:"chr2",4:"start2",5:"end2"},
                   inplace=True)
    assert domains.apply(lambda x: x["chr1"]==x["chr2"] and \
                                   x["start1"] == x["start2"] and \
                                   x["end1"] == x["end2"], axis="columns").all()
    assert (domains.end1 - domains.start1 > 0).all()
    domains.drop(["start2","chr2","end2"], axis="columns", inplace=True)
    domains.rename(columns={"chr1":"chr","start1":"start","end1":"end"},
                   inplace=True)
    domains["intervals"] = pd.arrays.IntervalArray.from_arrays(
        domains.start,domains.end,closed="both")
    domains.index = pd.MultiIndex.from_frame(domains[["chr", "intervals"]])
    return domains

# read compartments file as store as pd dataframe
def read_compartments(file):
    E1 = pd.read_csv(file,sep="\t",header=None)
    E1.rename(columns={0:"chr",1:"start",2:"end",3:"E1"}, inplace=True)
    E1.dropna(inplace=True)
    E1["intervals"] = pd.IntervalIndex.from_arrays(E1.start,E1.end,closed="both")
    E1.index=pd.MultiIndex.from_frame(E1[["chr","intervals"]])
    return E1

# old function, not really used
def statE1(domain,E1,f):
    overlap = E1.loc[domain.chr].index.overlaps(domain.intervals)
    e1 = E1.loc[domain.chr][overlap]
    if len(e1) == 0:
        print ("Warning, no E1 found for domain",domain)
        return 0
    return f(e1["E1"].values)

def plot_E1_from_domains_size_dependence():
    def getOverlapingDomainLength(e1, domains):
        overlap = domains.loc[e1.chr].index.overlaps(e1.intervals)
        domain = domains.loc[e1.chr][overlap]
        if len(domain) == 0:
            # print ("Warning, no domain found for e1 bin ",e1)
            return -1
        if len(domain) > 1:
            print("Warning, >1 domain found for e1 bin ", e1)
            raise
        l = domain.end.iloc[0] - domain.start.iloc[0]
        if l > 400000:
            return 400000 // 50000
        else:
            return l // 50000

    print("Getting domain length...")
    E1["domain_length"] = E1.apply(getOverlapingDomainLength, axis="columns", domains=domains)

    print("Plotting")
    ax = sns.violinplot(x="domain_length", y="E1", data=E1)
    ax.set(xlabel="Domain length, x50-kb")
    plt.axhline(y=0)
    plt.show()

def compartments_switch_at_domains_boundaries(domains, E1, useHash = False):
    def get_E1_near_boundaries(domain,E1,k,binsize):
        #get interval +/-k beens near boundary

        left_boundary_interval = pd.Interval(domain.start - k*binsize,
                                                   domain.start + (k+1)*binsize,
                                                   closed="neither")
        right_boundary_interval = pd.Interval(domain.end+1 - (k+1)*binsize,
                                                   domain.end+1 + k*binsize,
                                                   closed="neither")
        # print (E1.loc[domain.chr].index)
        # print (left_boundary_interval)
        overlap_left = E1.loc[domain.chr].index.overlaps(left_boundary_interval)
        overlap_left = E1.loc[domain.chr][overlap_left]

        overlap_right = E1.loc[domain.chr].index.overlaps(right_boundary_interval)
        overlap_right = E1.loc[domain.chr][overlap_right]

        # check overlap size
        if len(overlap_left) != (k*2 + 1):
            return None
        if len(overlap_right) != (k*2 + 1):
            return None

        return (overlap_left.E1.values.tolist() + overlap_right.E1.values.tolist())

        #get all e1 values vector
        #contactenate vectors for both boundaries
        #return

    k = 2 # how many bins near boundary to use
    hashfile = os.path.join("hashedData",
                            md5((domains_file+compartments_file+str(k)).encode()).hexdigest() + \
                            "."+compartments_switch_at_domains_boundaries.__name__+".dump")
    if useHash and os.path.exists(hashfile):
        domains = pickle.load(open(hashfile,"rb"))
    else:
        domains["E1_boundary"] = domains.apply(get_E1_near_boundaries,axis="columns",
                                           E1=E1,k=k,binsize=binsize)
        pickle.dump(domains,open(hashfile,"wb"))

    # E1_boundary = None when some of the E1 value missing
    # assert that this does not happen often
    assert pd.isna(domains["E1_boundary"]).sum() < (len(domains) / 10)
    print (domains.head())
    domains.query("end-start >= (@k+1)*@binsize-1",inplace=True)
    domains = domains[pd.notna(domains["E1_boundary"])]
    boundaries = domains["E1_boundary"].dropna().values.tolist()
    boundaries = np.array(boundaries)

    separate_boundaries = True
    # concat two boundaries
    if separate_boundaries:
        b1 = boundaries[:,:k*2+1]
        b2 = np.fliplr(boundaries[:,k*2+1:])
        print (boundaries)
        print ("---------")
        print (b1)
        print("---------")
        print (b2)
        print("---------")
        boundaries = np.vstack((b1,b2))
        print (boundaries)

    # set colormaps
    domainlengths = (domains.end-domains.start).values.tolist()
    norms = colors.BoundaryNorm(np.percentile(np.unique(domainlengths),
                                              range(0,100,20)),
                                ncolors=256)
    mapper = cm.ScalarMappable(norms,cmap="Greys")
    mapper.set_array(domainlengths)
    lencolors = list(map(mapper.to_rgba,domainlengths))
    if separate_boundaries:
        lencolors = lencolors*2

    E1_values = boundaries.flatten()
    E1mapper = colors.TwoSlopeNorm(vcenter=0,
                         vmin=np.percentile(E1_values[E1_values<0],5),
                         vmax=np.percentile(E1_values[E1_values>0],95))
    ax = sns.clustermap(boundaries,
                        row_cluster=True,
                        col_cluster=False,
                        metric="seuclidean",
                        #figsize=(max(1,(k*2+1)*2//10),
                        #         max(20,len(boundaries)//100)),
                        yticklabels=False,
                        row_colors=lencolors,
                        norm=E1mapper,
                        cmap="bwr"
                        )

    ax.ax_heatmap.axvline(x=k+0.5)
    if not separate_boundaries:
        ax.ax_heatmap.axvline(x=2*k+1+k+0.5)
        ax.ax_heatmap.axvline(x=2*k+1,ls="--")
    c = plt.colorbar(mapper, ax=ax.ax_col_dendrogram,
                     orientation="horizontal",
                     fraction=0.5)
    c.set_label("Domain size")

    plt.show()

# sie of the Hi-C bin. Should be same for E1 and domains
binsize = 25000

domains_file = "Domains/AcolNg_4.35.ann"
domains = read_domains(domains_file)
print (domains.head())
assert (domains.start % binsize == 0).all()

compartments_file = \
    "http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/eig/v3/AcolNg.v3.tpm.my.eig.bedGraph"
E1 = read_compartments(compartments_file)
print (E1.head())
assert (E1.start % binsize == 0).all()

print ("Starting analysis....")

# This func plots dependence of the E1-values of within-domain bins from domain size
# Motivated by the visual observation that long domains
# are predominantly located in the B-compartment

# plot_E1_from_domains_size_dependence()

# TODO: write a comment here
compartments_switch_at_domains_boundaries(domains, E1, useHash=True)