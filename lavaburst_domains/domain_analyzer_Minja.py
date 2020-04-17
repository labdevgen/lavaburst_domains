import pandas as pd
import seaborn as sns
from hashlib import md5
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import os
import pickle
import math
import numpy as np
import seaborn as sns
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

    c_res_mism_st = domains.start[domains.start % domains_resolution != 0].count()
    c_res_mism_end = domains.end[domains.end % domains_resolution != 0].count()
    if c_res_mism_st + c_res_mism_end != 0:
        print ("----WARNING! ",c_res_mism_st + c_res_mism_end," domain boundaries do not match resolution")
        print (c_res_mism_st," for starts, ",c_res_mism_end," for ends")
    domains["start"] = (domains["start"] // domains_resolution) * domains_resolution
    domains["end"] = (domains["end"] // domains_resolution) * domains_resolution

    assert (domains.start % domains_resolution == 0).all()
    assert (domains.end % domains_resolution == 0).all()

    domains["intervals"] = pd.arrays.IntervalArray.from_arrays(
        domains.start,domains.end,closed="both")
    domains.index = pd.MultiIndex.from_frame(domains[["chr", "intervals"]])
    domains["size"] = domains.end - domains.start
    assert np.all((domains["size"] > 0).values)
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

def statDomains(domains):
    sns.violinplot(x="chr", y="size", data=domains)
    plt.savefig("results/" + os.path.basename(domains_file) + ".sizedist.png")
    plt.clf()
    with open("results/" + os.path.basename(domains_file)+"_stats.txt","w") as fout:
        fout.write("Average domain size: " + str(domains["size"].median())+"\n")
        fout.write("Mean domain size: " + str(domains["size"].mean())+"\n")
        fout.write("Number of domains: " + str(len(domains))+ "\n")

# for each domain boundaries compute insulatory score
# insulatory score is simply an average obs/expected for rhombus with apex in domain boundary
# old version, use hicExplorer_scores for new domains
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

    #plt.savefig("results/"+os.path.basename(hic_file)+os.path.basename(domains_file)+".png")
    return domains

# for each domain boundaries read insulatory score from hicExplorer-based bedgraph
# old version, use hicExplorer_scores for new domains
def get_insulation_of_domain_boundaries_hicExplorer_scores(domains, scores_file, domains_resolution):
    def get_insulation(chr, pos, scores,resolution):
        if (chr,pos-resolution//2) in scores.index:
            return scores.loc[(chr,pos-resolution//2)]["score"]
        else:
            return np.nan

    scores = pd.read_csv(scores_file,sep="\t",header=None,names=["chr","start","end","score"])
    scores.set_index(keys=["chr","start"], inplace=True)
    assert scores.index.is_unique

    domains["insulation_l"] = domains.apply(lambda x: get_insulation(x.chr,x.start,scores,
                                                             domains_resolution),
                                            axis="columns")
    domains["insulation_r"] = domains.apply(lambda x: get_insulation(x.chr,x.end,scores,
                                                             domains_resolution),
                                            axis="columns")
    print (pd.isna(domains["insulation_l"]).sum()," left boundaries do not have defined insulation index")
    print(pd.isna(domains["insulation_r"]).sum(), " right boundaries do not have defined insulation index")
    domains.dropna(inplace=True)
    return domains

def plot_E1_from_domains_size_dependence(E1, length_bin = 25000, maxlength = 400000):
    def getOverlapingDomainLength(e1, domains, length_bin, maxlength):
        overlap = domains.loc[e1.chr].index.overlaps(e1.intervals)
        domain = domains.loc[e1.chr][overlap]
        if len(domain) == 0:
            # print ("Warning, no domain found for e1 bin ",e1)
            return -1
        if len(domain) > 1:
            #print("Warning, >1 domain found for e1 bin ", e1)
            return -1
        l = min(maxlength,domain["size"].iloc[0])
        return ((l // length_bin) * length_bin) // 1000 # return binned length in kb

    print("Getting domain length...")
    E1["domain_length"] = E1.apply(getOverlapingDomainLength, axis="columns",
                                   domains=domains, length_bin =length_bin, maxlength = maxlength)
    E1.query("domain_length != -1",inplace=True)
    print (E1["E1"].max(),E1["E1"].min(),E1["E1"].median())
    print("Plotting")
    ax = sns.violinplot(x="domain_length", y="E1", data=E1, palette="RdBu")
    ax.set_xlabel("Domain length, kb", fontsize=14)
    ax.set_ylabel("E1", fontsize=14)
    plt.axhline(y=0)
    plt.savefig("results/"+os.path.basename(domains_file)+".E1_vs_size.png",dpi=600)
    plt.clf()

def compartments_switch_at_domains_boundaries(domains, E1,
                                              score_file, domains_resolution,
                                              useHash = False):
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

    k = 3 # how many bins near boundary to use
    hashfile = os.path.join("hashedData",
                            md5((domains_file+compartments_file+str(k)).encode()).hexdigest() + \
                            "."+compartments_switch_at_domains_boundaries.__name__+".new.dump")
    if useHash and os.path.exists(hashfile):
        domains = pickle.load(open(hashfile,"rb"))
    else:
        domains.query("size >= 125000", inplace=True)
        domains = get_insulation_of_domain_boundaries_hicExplorer_scores(domains,scores_file=score_file,
                                                                         domains_resolution=domains_resolution)
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

    # whether to consider each boundary separately or draw whole domain
    # some plots work only for separate_boundaries = True
    separate_boundaries = True

    # concat two boundaries
    if separate_boundaries:
        b1 = boundaries[:,:k*2+1]
        b2 = np.fliplr(boundaries[:,k*2+1:])
        boundaries = np.vstack((b1,b2))
        boundary_index = boundaries.shape[1] // 2 + 1

    ######################## dendrogramm #################
    print("Drawing boundaries dendrogramm")

    # set colormaps
    domainlengths = (domains.end-domains.start).values.tolist()
    domain_insulations = np.hstack((domains.insulation_l.values,domains.insulation_r.values))

    length_norms = colors.BoundaryNorm(np.percentile(np.unique(domainlengths),
                                              range(0,100,20)),
                                ncolors=256)
    insulation_norms = colors.BoundaryNorm(np.percentile(np.unique(domain_insulations),
                                              range(0,100,10)),
                                ncolors=256)

    length2color_mapper = cm.ScalarMappable(length_norms,cmap="Greys")
    length2color_mapper.set_array(domain_insulations)

    insulation2color_mapper = cm.ScalarMappable(insulation_norms,cmap="RdBu")

    inscolors_l = list(map(insulation2color_mapper.to_rgba,domains.insulation_l))
    inscolors_r = list(map(insulation2color_mapper.to_rgba,domains.insulation_r))

    lencolors = list(map(length2color_mapper.to_rgba,domainlengths))

    if separate_boundaries:
        colors_list = [lencolors*2,inscolors_l+inscolors_r]
    else:
        colors_list = [lencolors,inscolors_l,inscolors_r]

    E1_values = boundaries.flatten()
    E1mapper = colors.TwoSlopeNorm(vcenter=0,
                         vmin=np.percentile(E1_values[E1_values<0],5),
                         vmax=np.percentile(E1_values[E1_values>0],95))

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage


    if separate_boundaries:
        averaged = np.vstack((np.average(boundaries[:,:boundary_index],axis=1),
                              np.average(boundaries[:,boundary_index+1:],axis=1))).T
        assert len(averaged)==len(boundaries)
        """
        distances = pdist(averaged, metric=lambda x,y: abs((x[0]-x[1])-(y[0]-y[1])))
        # distances = np.nan_to_num(distances,nan=0)
        link = linkage(distances, optimal_ordering=True)
    else:
        link = linkage(pdist(boundaries), optimal_ordering=True)

    ax = sns.clustermap(boundaries,
                        row_cluster=True,
                        row_linkage=link,
                        col_cluster=False,
                        #metric="seuclidean",
                        #figsize=(max(1,(k*2+1)*2//10),
                        #         max(20,len(boundaries)//100)),
                        yticklabels=False,
                        row_colors=colors_list,
                        norm=E1mapper,
                        cmap="bwr"
                        )
    ax.ax_heatmap.axvline(x=k+0.5)
    if not separate_boundaries:
        ax.ax_heatmap.axvline(x=2*k+1+k+0.5)
        ax.ax_heatmap.axvline(x=2*k+1,ls="--")
    c = plt.colorbar(length2color_mapper, ax=ax.ax_col_dendrogram,
                     orientation="horizontal",
                     fraction=0.5)
    c.set_label("Domain size")

    plt.show()
    """

    plt.clf()

    # next plots are only for separate boundaries
    if not separate_boundaries:
        return

    ################## scatter plot #####################
    print("Drawing E1-diff vs insulation scatterplot")
    # on X-axes E1 difference
    # on Y-axes insulation score
    sns.scatterplot(x=np.subtract(averaged[:,0],averaged[:,1]),
                    y=domains.insulation_l.values.tolist()+domains.insulation_r.values.tolist()
                )
    #plt.show()

    ################# examples of insulated domains w/o E1 change #############
    abs_E1diff = np.abs(np.subtract(averaged[:,0],averaged[:,1]))
    E1diff_threashold = np.percentile(abs_E1diff,25)
    insulation_threashold = np.percentile(domains.insulation_l.values.tolist()+\
                                      domains.insulation_r.values.tolist(),
                                      25)
    assert len(abs_E1diff) % 2==0
    mask_l = (domains.insulation_l.values < insulation_threashold) & \
         (abs_E1diff[:len(abs_E1diff)//2] < E1diff_threashold)
    mask_r = (domains.insulation_r.values < insulation_threashold) & \
             (abs_E1diff[len(abs_E1diff)//2:] < E1diff_threashold)

    domain_colors = np.array(["255,255,255"]*len(domains))
    domain_colors[mask_l] = "128,128,0" # Olive
    domain_colors[mask_r] = "0,0,255" # Blue
    domain_colors[mask_r & mask_l] = "255,255,0" # Yellow
    domains["color"] = domain_colors

    examples = domains[mask_l | mask_r]
    def writeJuicerAnnotation(d,f,color):
        f.write("\t".join(map(str,[d.chr,d.start,d.end,d.chr,d.start,d.end,color]))+"\n")

    basename ="results/"+os.path.basename(domains_file)
    for mask, file, color in zip([np.logical_not(mask_l) & mask_r,
                                  np.logical_not(mask_r) & mask_l,
                                    (mask_l & mask_r)],
                                ["r","l","b"],
                                ["0,0,255","128,128,0","255,255,0"]):
      with open(basename+"insulated.examples"+file+".ann","w") as fout:
        domains.iloc[mask,:].apply(writeJuicerAnnotation, axis="columns", f=fout, color=color)

    raise
    pd.concat(
        (pd.DataFrame([domains.chr,
                      domains.start,
                      domains.start+binsize,
                      domains.insulation_l]).rename(columns={"insulation_l":"insulation"}),
        pd.DataFrame([domains.chr,
                      domains.end+1,
                      domains.end+1+binsize,
                      domains.insulation_r]).rename(columns={"insulation_r":"insulation"})
         ), ignore_index=True).to_csv(basename+"ins.bedGraph",sep="\t",index=False)

# sie of the Hi-C bin. Should be same for E1 and domains
domains_resolution = 5000

# hic_file = "hics/AcolNg_V3.hic.25000.oe.1000000.MB"
# domains_file = "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hicExplorer_TADs/results/AalbS2_V4.1000.hic_5000.h5.delt.0.05_domains.2D"
domains_file = "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hicExplorer_TADs/results/AatrE3_V4.1000.hic_5000.h5.delt.0.05_domains.2D"
scores_file = domains_file.replace("_domains.2D","_score.bedgraph")

compartments_file = \
    "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/eig/v3/AatrE3_V3.tpm.my.eig.bedGraph"

print (domains_file)
domains = read_domains(domains_file)

# plot violin plot of domain sizes, print mean and average to file
# statDomains(domains)

E1 = read_compartments(compartments_file)
# assert (E1.start % binsize == 0).all()

print ("Starting analysis....")

# This func plots dependence of the E1-values of within-domain bins from domain size
# Motivated by the visual observation that long domains
# are predominantly located in the B-compartment

#plot_E1_from_domains_size_dependence(E1)


# This function will add insulation score for each genomic boundary
# technically this will add insulation_r / l fields to domain dframe
get_insulation_of_domain_boundaries_hicExplorer_scores(domains,scores_file=scores_file,
                                                       domains_resolution=domains_resolution)
raise

for offset in [1,2]:
   for dist in range(offset+1,7):
       print (offset,dist)
       get_insulation_of_domain_boundaries(domains, offset=offset, dist=dist)
get_insulation_of_domain_boundaries(domains)

compartments_switch_at_domains_boundaries(domains, E1, useHash=True)