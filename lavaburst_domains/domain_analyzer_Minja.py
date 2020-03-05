import pandas as pd
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
    return domains

def read_compartments(file):
    E1 = pd.read_csv(file,sep="\t",header=None)
    E1.rename(columns={0:"chr",1:"start",2:"end",3:"E1"}, inplace=True)
    E1["intervals"] = pd.IntervalIndex.from_arrays(E1.start,E1.end,closed="both")
    E1.index=pd.MultiIndex.from_frame(E1[["chr","intervals"]])
    print (E1.head())
    print (E1.loc["X"])
    return E1

def statE1(domain,E1,f):
    # SUPER SLOW =)
    e1 = E1[E1.loc[domain.chr].intervals.overlaps(domain.intervals)]
    if len(e1) == 0:
        print ("Warning, no E1 found for domain",domain)
        return 0
    return f(e1["E1"].values)

domains = read_domains("Domains/AcolNg_4.35.ann")
print (domains.head())

E1 = read_compartments(
    "http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/eig/v3/AcolNg.v3.tpm.my.eig.bedGraph")
print (E1.head())

domains["averageE1"] = domains.apply(statE1,axis="columns",E1=E1,f=np.median)

colors = {"X":"r","2L":"g","2R":"b","3R":"k","3L":"k"}
plt.scatter(domains.end - domains.start, domains.averageE1,
            color = domains["chr"].apply(lambda x:colors [x]))
plt.show()