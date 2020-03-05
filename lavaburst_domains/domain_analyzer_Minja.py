import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

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

def read_compartments(file):
    E1 = pd.read_csv(file,sep="\t",header=None)
    E1.rename(columns={0:"chr",1:"start",2:"end",3:"E1"}, inplace=True)
    E1["intervals"] = pd.IntervalIndex.from_arrays(E1.start,E1.end,closed="both")
    E1.index=pd.MultiIndex.from_frame(E1[["chr","intervals"]])
    return E1

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

def compartments_switch_at_domains_boundaries():
    def get_E1_near_boundaries(domain,E1,k):
        #get interval +/-k beens near boundary
        #get all e1 values vector
        #contactenate vectors for both boundaries
        #return
        pass

    domains["E1_boundary"] = domains.apply(get_E1_near_boundaries,axis="columns",
                                           E1=E1,k=3)


binsize = 25000
domains_file = "Domains/AcolNg_4.35.ann"
domains = read_domains(domains_file)
print (domains.head())

E1 = read_compartments(
    "http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/eig/v3/AcolNg.v3.tpm.my.eig.bedGraph")
print (E1.head())

plot_E1_from_domains_size_dependence()