import pandas as pd
import seaborn as sns
from hashlib import md5
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import os
import pickle
import math
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

# for each domain boundaries compute insulatory score
# insulatory score is simply an average obs/expected for rhombus with apex in domain boundary

def get_insulation_of_domain_boundaries(domains, offset = 1, dist = 3):
    def get_insulation(chr, position, binsize, hic_data):
        if position - (dist+1)*binsize <= 0:
            return np.nan
        oes = []
        for start in range(position - dist*binsize,
                           position - offset*binsize,
                           binsize):
            for end in range(position + offset*binsize,
                             position + dist*binsize,
                             binsize):
                oe = np.nan
                try:
                    oe = hic_data.loc[chr,start,end].oe
                except KeyError:
                    pass
                oes.append(oe)
        if sum(np.isnan(oes)) == len(oes):
            return np.nan
        return np.nanmean(oes)

    assert dist > offset
    hic_data = pd.read_csv(hic_file,sep="\t",header=None).rename(
                                    columns={0:"chr",1:"start",2:"end",3:"oe"})
    hic_data.oe = hic_data.oe.apply(math.log)
    hic_data.index = pd.MultiIndex.from_frame(hic_data[["chr", "start","end"]])
    #hic_data.drop(columns=["chr","start","end"],inplace=True)
    hic_data.dropna(axis=0, inplace=True)

    domains["insulation_l"] = domains.apply(lambda x: get_insulation(x.chr,x.start,
                                                             binsize,hic_data),
                                     axis="columns")
    domains["insulation_r"] = domains.apply(lambda x: get_insulation(x.chr,x.end + 1,
                                                             binsize,hic_data),
                                            axis="columns")
    domains.dropna(inplace=True)

    sample = hic_data.query( "(end - start <= 2*@dist*@binsize) & " +
                             "(end - start >= 2*@offset*@binsize)")
    source = ["TAD left boundary"]*len(domains)+\
             ["TAD right boundary"]*len(domains)+\
             ["Genome average"]*len(sample)
    insulation = np.hstack((domains["insulation_l"].values,
                                 domains["insulation_r"].values,
                                 sample.oe.values)).flatten()
    plot_data = pd.DataFrame({"region":source, "log(obs/exp)":insulation})
    ax = sns.boxplot(x="region", y="log(obs/exp)",data=plot_data)
    ax.set(ylabel="log(obs/exp)")
    plt.axhline(0)
    print (np.median(domains["insulation_l"].values),
                            np.median(sample.oe.values),
           np.median(domains["insulation_l"].values) - np.median(sample.oe.values)
           )
    #plt.show()

    plt.savefig("results/"+os.path.basename(hic_file)+os.path.basename(domains_file)+".png")
    return domains

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
    domains = get_insulation_of_domain_boundaries(domains,E1)
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
        boundaries = np.vstack((b1,b2))

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

hic_file = "hics/AcolNg_V3.hic.25000.oe.1000000.MB"
domains_file = "Domains/AcolNg_4.35.ann"
compartments_file = \
    "http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/eig/v3/AcolNg.v3.tpm.my.eig.bedGraph"

domains = read_domains(domains_file)
print (domains.head())
assert (domains.start % binsize == 0).all()

E1 = read_compartments(compartments_file)
print (E1.head())
assert (E1.start % binsize == 0).all()

print ("Starting analysis....")

# This func plots dependence of the E1-values of within-domain bins from domain size
# Motivated by the visual observation that long domains
# are predominantly located in the B-compartment

# plot_E1_from_domains_size_dependence()

# This function will add insulation score for each genomic boundary
# technically this will add insulation_r / l fields to domain dframe

#for offset in [1,2]:
#    for dist in range(offset+1,7):
#        print (offset,dist)
#        get_insulation_of_domain_boundaries(domains, offset=offset, dist=dist)
get_insulation_of_domain_boundaries(domains)

# compartments_switch_at_domains_boundaries(domains, E1, useHash=True)