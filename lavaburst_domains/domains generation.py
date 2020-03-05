#!/usr/bin/env python
# coding: utf-8

# In[1]:


from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import lavaburst
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import makedirs
import time
import random
from itertools import product
from scipy.signal import argrelmax,argrelmin

binsize = 25000


# In[2]:



"""
In order to get matrices, I used this script:


JUICER=/home/al/Juice/juicer.jar
StringVal="X 2R 2L 3R 3L"
for i in $StringVal
do
    echo $i
    java -jar $JUICER dump observed -d KR http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg_V3.1000.hic $i $i BP 25000 /home/al/Vp/Matrix25000/AcolNg1000_${i}_25000.txt
done
"""


# In[3]:


species = ['Aalb','Aatr','AcolNg','Amer','Aqmac','Aste','AmerA']
chroms = ['X','2R','2L','3R','3L']
# species = ['AcolNg']
# chroms = ['X']

ins_track = []
specific_minima = []
im_neighbor = []
im_exact_intersect = []

for sp, chrom in product(species,chroms):
    A = np.loadtxt('/home/al/Vp/Matrix25000/'+sp+'1000_' + chrom + '_25000.txt')
    A[np.isnan(A)] = 0
    l = []
    #Only four isolation window values are allowed in this code.
    for i in [2,3,4,5]:
        ins = lavaburst.utils.insul_diamond(A,i)
        ins = pd.DataFrame(ins)
        ins[1] = ins.index
        ins[1] *= binsize
        ins[2] = ins[1] + binsize -1
        ins = ins[[1,2,0]]
        ins.columns = [0,1,2]
        ins.insert(0, 0, chrom, allow_duplicates=True)
        ins.columns = [0,1,2,3]
        fullins = np.array(ins)
        ins_track.append([sp,i,fullins])
        
        ins[3] = ins[ins.columns[3]][(ins[ins.columns[3]].shift(1) > 0.02 + ins[ins.columns[3]]) & (ins[ins.columns[3]].shift(-1) > 0.02 + ins[ins.columns[3]])]
        ins[3].loc[ins[3].notna()] = 1
        ins.fillna(0,inplace=True)
        ins = ins.reset_index(drop=True)
        
        part_mins = ins.loc[ins[3]==1]
        part_mins = part_mins[[0,1,2]].reset_index(drop=True)
#         part_mins = np.array(part_mins)
        specific_minima.append([sp,i,part_mins])
    
        l.append(ins[3])
    ins = pd.concat(l,axis=1)
    ins.columns = [i for i in range(0,ins.shape[1])]
    ins[len(ins.columns)] = ins.sum(axis=1)
#     print(ins)
    a = np.array(ins[4])
    for i in range(1,len(a)):
        if a[i]==3:
#             print(a[i-1],a[i])
            if a[i-1]==1:
                a[i] += 1
                a[i-1] -= 1
        elif a[i]==2:
            if a[i-1]==2:
                a[i] += 2
                a[i-1] -= 2 
        elif a[i]==1:
            if a[i-1]==3:
                a[i] -= 1
                a[i-1] += 1
    ins1 = ins.loc[ins[4]>=4]
    ins[4] = a
    ins = ins.loc[ins[4]>=4]
    
    ins[0] = ins.index
    ins1[0] = ins1.index
    ins[0] = ins[0]*binsize
    ins1[0] = ins1[0]*binsize
    ins[1] = ins[0] + binsize -1
    ins1[1] = ins1[0] + binsize -1
    ins.insert(0, 0, chrom, allow_duplicates=True)
    ins.columns = [i for i in range(0,ins.shape[1])]
    ins1.insert(0, 0, chrom, allow_duplicates=True)
    ins1.columns = [i for i in range(0,ins1.shape[1])]
    ins = ins[[0,1,2]]
    ins1 = ins1[[0,1,2]]
    ins = np.array(ins)
    ins1 = np.array(ins1)
    im_neighbor.append([sp,ins])
    im_exact_intersect.append([sp,ins1])

ins_track = np.array(ins_track)
specific_minima = np.array(specific_minima)
im_neighbor = np.array(im_neighbor)
im_exact_intersect = np.array(im_exact_intersect)

ws = [2,3,4,5]
"""
Generates an insulation track for the specified window sizes.
"""
for sp,w in product(species,ws):
    fi2 = ins_track[(ins_track[:,0]==str(sp))&(ins_track[:,1]==w)]
    fi2 = fi2[:,2]
    fi3 = np.concatenate(fi2,axis=0)
    fi3 = pd.DataFrame(fi3)
    fi3 = fi3.sort_values([0,1,2]).reset_index(drop=True)
#     np.savetxt('/home/al/Vp/Ins/'+sp+'_ins'+str(w)+'.bedGraph', fi3, delimiter='\t', newline='\n',fmt='%s')

"""
"Extended" intersection of insulation minima between four window sizes.

For example, if for window sizes w = [2,4,5], 
insulation minimums are in the bin with the coordinate "x", 
but for window size w = [3], minimums are in the bin with the coordinate "x-1", 
such a missmatch is not discarded.

The bin "x" is considered to have a minimum that intersects for all four windows.
The scheme of the algorithm:
(1|3 = 0|4)
(2|2 = 0|4)
(3|1 = 4|0)

"""
for sp in species:
    immw = im_neighbor[im_neighbor[:,0]==str(sp)]
    immw = immw[:,1]
    immw2 = np.concatenate(immw,axis=0)
#     np.savetxt('/home/al/Vp/IM/'+sp+'_IM.bed', immw2, delimiter='\t', newline='\n',fmt='%s')

#Exact intersection of insulation minima between several window sizes.
# for sp in species:
#     imm = im_exact_intersect[im_exact_intersect[:,0]==str(sp)]
#     imm = imm[:,1]
#     imm2 = np.concatenate(imm,axis=0)
# #     np.savetxt('/home/al/Vp/IMT/'+sp+'_IMwout.bed', imm2, delimiter='\t', newline='\n',fmt='%s')
  
# #Create insulation minima for a specific window size.
# for sp,w in product(species,ws):
#     part_m1 = specific_minima[(specific_minima[:,0]==str(sp))&(specific_minima[:,1]==w)]
#     part_m1 = part_m1[:,2]
#     part_m2 = np.concatenate(part_m1,axis=0)
#     part_m2 = pd.DataFrame(part_m2)
#     part_m2 = part_m2.sort_values([0,1,2]).reset_index(drop=True)

# #     np.savetxt('/home/al/Vp/IMT29/'+sp+'_IM'+str(w)+'.bed', part_m2, delimiter='\t', newline='\n',fmt='%s')


# In[4]:


species = ['Aalb','Aatr','AcolNg','Amer','Aqmac','Aste','AmerA']
chroms = ['X','2R','2L','3R','3L']
species = ['Aste']
chroms = ['X']
sp_delta = []
"""
Generates a delta vector according to:
https://www.ncbi.nlm.nih.gov/pubmed/26030525

Insulation diamond window: 5 bins.
Delta window: 5 bins.

Delta vector = (mean (2 left bins)) - (mean (2 right bins)).

"""
for sp, chrom in product(species,chroms):
    ins = pd.read_csv('/home/al/Vp/Ins/'+sp+'_ins5.bedGraph',sep='\t',header = None)
    f_vals = ins.loc[(ins[0] == chrom)].reset_index(drop=True)
    f_vals1 = pd.DataFrame(f_vals[3])
    f_vals1[1] = f_vals.index
    f_vals1 = f_vals1[[1,3]]
    f_vals1.columns = [0,1]
    
    insmins = f_vals1.copy()
    insmins[1] = f_vals1[f_vals1.columns[1]][(f_vals1[f_vals1.columns[1]].shift(1) > f_vals1[f_vals1.columns[1]]) 
                                          & (f_vals1[f_vals1.columns[1]].shift(-1) > f_vals1[f_vals1.columns[1]])]
    insmins = insmins.dropna()
    insmins = insmins.reset_index(drop=True)
    insmins = insmins[0].values

    f_vals1[2] = f_vals1[1].shift(1)
    f_vals1[3] = f_vals1[1].shift(2)
    f_vals1[4] = f_vals1[1].shift(3)
    f_vals1[5] = f_vals1[1].shift(-1)
    f_vals1[6] = f_vals1[1].shift(-2)
    f_vals1[7] = f_vals1[1].shift(-3)
    f_vals1[8] = f_vals1[[2,3]].mean(axis=1,skipna = False)
    f_vals1[9] = f_vals1[[5,6]].mean(axis=1,skipna = False)

    f_vals1 = f_vals1[[0,1,8,9]]
    f_vals1.columns = [0,1,2,3]
    f_vals1[4] = f_vals1[2] - f_vals1[3]
    f_vals1.drop(f_vals1.head(3).index,inplace=True)
    f_vals1.drop(f_vals1.tail(3).index,inplace=True)
    f_vals1[5] = f_vals1[4].shift(-1)
    f_vals1[6] = 0
    f_vals1[6].loc[f_vals1.loc[f_vals1[4]>=0].loc[f_vals1[5]<=0].index] += 1
    f_vals2 = f_vals1.loc[f_vals1[6]==1]
    mins = f_vals2[0]
    startend = pd.DataFrame(data = {0:[0,max(f_vals1[0])]})
    startend.index = startend[0].tolist()
    mins = mins.append(startend).sort_values([0])
    mins[1] = mins[0].shift(-1)
    mins[2] = mins[0].shift(-2)
    mins.drop(mins.tail(2).index,inplace=True)
    
    armin = np.array(mins)
    arref = np.array(f_vals1)
    ar2 = []
    for m in armin:
        m0 = m[0]
        m1 = m[1]
        m2 = m[2]
        leftmax = max(arref[(arref[:,0]>=m0) & (arref[:,0]<=m1)][:,4])
        if (leftmax < 0.05):
            continue
        rightmin = min(arref[(arref[:,0]>=m1) & (arref[:,0]<=m2)][:,4])
        if (rightmin > -0.05):
            continue
        #boundary_strenght
        bs = leftmax-rightmin
        bs_sn = 1. / (1 + np.exp(-bs)) * 2 - 1
        if (bs_sn < 0.1):
            continue
        if m1+1 in insmins:
            m1 += 1
        ar2.append([m1,bs,bs_sn])
    ar2 = pd.DataFrame(ar2)

    ar2.loc[:,0] *= binsize
    ar2[3] = ar2[0]+binsize-1
    ar3 = ar2[[0,3,2]]
    ar3.insert(0, 0, chrom, allow_duplicates=True)
    ar3.columns = [0,1,2,3]
    ar3 = ar3.astype({0:str,1:int,2:int})
    sp_delta.append([sp,ar3])
sp_delta = np.array(sp_delta)

for sp in species:
    sp_delta1 = sp_delta[sp_delta[:,0]==str(sp)]
    sp_delta2 = sp_delta1[:,1]
    sp_delta3 = np.concatenate(sp_delta2,axis=0)
#     np.savetxt('/home/al/Vp/Delta/'+sp+'_delta.bed', sp_delta3, delimiter='\t', newline='\n',fmt='%s')


# In[5]:


"""
Requirements:
- lavaburst package https://github.com/nvictus/lavaburst ;
- Dense HiC matrix;
- File with restricted starts/ends for each chromosome;
- Insulation vector;
- Insulation minima (intersected* between four window sizes);
- Delta vector.
"""



def consensus_domains(segments, weights):
    """
    Returns consensus list of nonoverlapping segments.
    Segments are 2-tuples given as half-open intervals [a,b).

    """
    occ = defaultdict(int)
    for d, w in zip(segments, weights):
        occ[d] += w
    # map each domain to its closest non-overlapping predecessor
    M = len(segments)
    seg = pd.DataFrame(segments)
    prev = np.zeros(M, dtype=int)
    for i in range(M-1, -1, -1):
        d = segments[i]
        j = i - 1
        while j > -1:
            if segments[j][1] <= d[0]: 
                prev[i] = j
                break
            j -= 1
    pr = pd.DataFrame(prev)
    # weighted interval scheduling dynamic program
    score = np.zeros(M, dtype=float)
    for i in range(1, M):
        d = segments[i]
        s_choose = score[prev[i]] + occ[d]
        s_ignore = score[i-1]
        score[i] = max(s_choose, s_ignore)
    consensus = []
    j = M - 1
    while j > 0:
        if score[j] != score[j-1]:
            consensus.append(segments[j])
            j = prev[j]
        else:
            j -= 1
    return consensus[::-1], max(score)


def sets_generator (A, g1, g2, gstep):
    """
    g1 - min gamma, g2 - max gamma.
    """
    seglist = []
    for g in np.arange(g1,g2,gstep): 
        S1 = lavaburst.scoring.modularity_score(A, gamma=g)
        model1 = lavaburst.SegModel(S1)
        segs,sc = model1.optimal_segmentation(return_score=True)
        sc = np.array(sc,ndmin=2).transpose()
        segs = np.concatenate([segs,sc],axis=1)
        seglist.append(segs)
    return seglist




def column_insul(col):
    
    """
    The function checks if the boundaries of the TAD intersect with
    minima of insulation or delta minima for each TAD border.

    It returns the insulation values at the border and 
    the intersection scores with the corresponding minima.
    """
    
    insul_minima_score = 0.25
    delta_minima_score = 0.5
    
    df = pd.DataFrame()
    df['ins_values'] = f_vals1.loc[col][1].reset_index(drop=True)
    df[2] = 0
    df.fillna(0, inplace = True)
    df = df.set_index(col.values,drop=True)


    df['ins_minima'] = 0
    bm = np.intersect1d(col,insmins[:,0])
    df.loc[bm,'ins_minima'] = insul_minima_score
    


    df['delta_minima'] = 0
    bd = np.intersect1d(col,deltas[:,0])
    df.loc[bd,'delta_minima'] = delta_minima_score
    
#     #Don't use.
#     am = np.intersect1d(col-1,insmins[:,0])
#     cm = np.intersect1d(col+1,insmins[:,0])
#     df.loc[(am+1),3] = 0.25
#     df.loc[(cm-1),3] = 0.25
#     ad = np.intersect1d(col-1,deltas[:,0])
#     cd = np.intersect1d(col+1,deltas[:,0])
#     df.loc[(ad+1),4] = 2
#     df.loc[(cd-1),4] = 2
    return df['ins_values'].values, df['delta_minima'].values, df['ins_minima'].values


def modul_score(SA,table):
    """
    The function assigns scores to the domain 
    if the domain belongs to the corresponding percentile interval. 

    The function takes an array of segment aggregation (modularity) score,
    a list of lists [[percentile start value, percentile end value, score]]
    and returns corresponding modularity TAD score. 

    i[0],i[1] - left, right percentile boundary values of modularity;
    i[2] - percentile interval score.
    """
    l = []
    for i in table:
        sc = np.copy(SA)
        sc[(sc>=i[0])&(sc<i[1])] = i[2]
        sc[(sc!=i[2])]=0
        l.append(sc)
    sc1 = sum(l)
    return sc1

species = ['Aalb', 'Aatr', 'AcolNg', 'Amer', 'Aqmac', 'AmerA']
chroms = ['X', '2R', '2L', '3R', '3L']
species = ['Aste']
chroms = ['3L']

"""
Set the range of gamma parameter.
"""
gg1 = [4]
gg2 = [x for x in range(15, 35, 3)]
gg2 = [x for x in range(35, 90, 5)]
gg2 = gg2 + [x for x in range(80, 180, 20)]
gg2 = gg2 + [x for x in range(180, 240, 30)]
gg2 = gg2 + [x for x in range(240, 600, 100)]
gg2 = [105] 
ggstep = [1]

"""
Set the number of domain sets with different average tad lengths.
"""
sizes = [7.2, 5.6, 4.35]

meta = []
start_time = time.process_time()
for sp, chrom in product(species,chroms):
    #Restricted start/end for each chromosome.
    rchr = pd.read_csv('/home/al/Vp/RSizes/'+sp+'.sizes',sep='\t',header = None)
    rchr1 = rchr.loc[(rchr[0] == chrom)]
    rchr1 = rchr1.reset_index(drop=True)
    r_start = int(rchr1[1][0]/binsize)
    r_end = int(rchr1[2][0]/binsize)
    assert r_start > 0

    ins = pd.read_csv('/home/al/Vp/Ins/'+sp+'_ins4.bedGraph',sep='\t',header = None)
    mins = pd.read_csv('/home/al/Vp/IM/'+sp+'_IM.bed',sep='\t',header = None) 
    deltas = pd.read_csv('/home/al/Vp/Delta/'+sp+'_delta.bed',sep='\t',header = None) 

    B = pd.read_csv('/home/al/Vp/Matrix25000/'+sp+'1000_'+chrom+'_25000.txt',delimiter='\t',header = None)
    B.drop(B.columns[len(B.columns)-1], axis=1, inplace=True)
    assert B.shape[0] == B.shape[1]
    A = np.array(B)    
    A[np.isnan(A)] = 0

    f_vals = ins.loc[(ins[0] == chrom)].reset_index(drop=True)
    f_vals1 = pd.DataFrame(f_vals[3])
    f_vals1[1] = f_vals.index
    f_vals1 = f_vals1[[1,3]]
    f_vals1.columns = [0,1]
    ff = np.array(f_vals1)

    nmins = f_vals1.copy()
    nmins[1] = f_vals1[f_vals1.columns[1]][(f_vals1[f_vals1.columns[1]].shift(1) > f_vals1[f_vals1.columns[1]]) 
                                          & (f_vals1[f_vals1.columns[1]].shift(-1) > f_vals1[f_vals1.columns[1]])]
    nmins = nmins.dropna()
    nmins = nmins.reset_index(drop=True)

    insmins = mins.loc[(mins[0] == chrom)]
    insmins = pd.DataFrame(insmins[1])
    insmins.columns = [0]
    insmins /= binsize
    insmins = f_vals1.loc[insmins[0].values].reset_index(drop=True)
    insmins = np.array(insmins)

    deltas = deltas.loc[(deltas[0] == chrom)]
    deltas = pd.DataFrame(deltas[[1,3]]).reset_index(drop=True)
    deltas.columns = [0,1]
    deltas[0] /= binsize
    deltas = np.array(deltas)
    
    """
    Alternative domains between two insulation minima.
    Domain quality is assessed 
    by the value of a strong maximum of insulation inside 
    minus the average value of insulation at the borders
    and by the absence of delta minima inside.
    Experimental feature.
    """
    
    dmins = nmins.copy()
    dmins[1] = dmins[0].shift(-1)
    dmins = dmins.dropna()
    dmins[2] = 0
    dmins = dmins.loc[dmins[0]>=r_start].reset_index(drop=True)
    dmins = dmins.loc[dmins[1]<r_end]
    dmins = dmins.astype({0:int,1:int})
    dmins1 = dmins[[0,1]].apply(column_insul,axis=0)
    dmins1 = np.transpose(np.concatenate(dmins1))
    dmins1 = pd.DataFrame(dmins1)
    dmins1[6] = dmins1[[2,5]].sum(axis=1)
    dmins1[7] = dmins1[[1,4]].sum(axis=1)
    dmins1 = pd.concat([dmins[[0,1,2]],dmins1],axis=1)
    dmins1.columns = [x for x in range(0,dmins1.shape[1])]
    dmins1 = dmins1[[0,1,2,3,6,9,10]]
    dmins1.columns = [x for x in range(0,dmins1.shape[1])]

    aaar11 = np.array(dmins1)
    aaar22 = []
    for i in aaar11:
        x = i[0]
        y = i[1]
        ix = i[3]
        iy = i[4]
        ibmean = (ix+iy)/2

        try:
            ff1 = ff[(ff[:,0]>=x) & (ff[:,0]<=y)]
            maxins1 = max(ff1[argrelmax(ff1[:,1])][:,1])
            maxins1 = maxins1 - ibmean
        except IndexError:
            maxins1 = 0
        except ValueError:
            maxins1 = 0
            
        try:
            dlt_inside = max(deltas[(deltas[:,0]>x+1) & (deltas[:,0]<y-1)][:,1])
        except IndexError:
            dlt_inside = 0
        except ValueError:
            dlt_inside = 0

        if dlt_inside > 0.15:
            d_insd = -9
        elif (dlt_inside >0) & (dlt_inside<=0.15):
            d_insd = -4
        else:
            d_insd = 0
        i = np.append(i,[0,d_insd,maxins1]).tolist()
        aaar22.append(i)

    dmins = np.array(aaar22)
    small = dmins[(dmins[:,1]-dmins[:,0])<4]
    small75p = np.percentile(small[:,9],75)
    small90p = np.percentile(small[:,9],90)
    small75 = small[(small[:,9]>small75p)&(small[:,9]<=small90p)]
    small90 = small[small[:,9]>small90p]
    small75[:,9] = 0.5
    small90[:,9] = 2

    big = dmins[(dmins[:,1]-dmins[:,0])>=4]
    big90p = np.percentile(big[:,9],90)
    big80p = np.percentile(big[:,9],80)
    big70p = np.percentile(big[:,9],70)
    big60p = np.percentile(big[:,9],60)
    big50p = np.percentile(big[:,9],50)
    big90 = big[big[:,9]>big90p]
    big80 = big[(big[:,9]>big80p)&(big[:,9]<=big90p)]
    big70 = big[(big[:,9]>big70p)&(big[:,9]<=big80p)]
    big60 = big[(big[:,9]>big60p)&(big[:,9]<=big70p)]
    big50 = big[(big[:,9]>big50p)&(big[:,9]<=big60p)]
    big50[:,9] = 0.5
    big60[:,9] = 0.75
    big70[:,9] = 1.5
    big80[:,9] = 2.5
    big90[:,9] = 3

    dmins = np.concatenate([small75,small90,big50,big60,big70,big80,big90],axis=0)

    dmins = pd.DataFrame(dmins)
    
    dmins[10] = 0
    dmins[11] = dmins[[5,6,7,8,9,10]].sum(axis=1)
    dmins[11].iloc[dmins.loc[dmins[11] == 0].index]=0.001
    dmins = dmins[[0,1,2,3,4,5,6,7,8,10,9,11]]
    dmins.columns = range(dmins.shape[1])
    
    for g1,g2,gstep in product(gg1,gg2,ggstep):
        seglist= sets_generator(A,g1,g2+0.01,gstep)

        updsegs = []
        for seg in seglist:
            seg = seg[seg[:,0]>=r_start]
            seg = seg[seg[:,1]<r_end]
            seg = seg[(seg[:,1]-seg[:,0])>=2]
            try:
                start = seg[:,0].min()
            except ValueError:
                continue
            end = seg[:,1].max()
            seg = pd.DataFrame(seg)
            seg = seg.astype({0:int,1:int})

            df1 = seg[[0,1]].apply(column_insul,axis=0)
            df1 = np.transpose(np.concatenate(df1))
            df1 = pd.DataFrame(df1)
            
            df1[6] = df1[[2,5]].sum(axis=1)
            df1[7] = df1[[1,4]].sum(axis=1)
            
            df1 = pd.concat([seg[[0,1,2]],df1],axis=1)
            df1.columns = [x for x in range(0,df1.shape[1])]
            df1 = df1[[0,1,2,3,6,9,10]]
            df1.columns = [x for x in range(0,df1.shape[1])]
            ar1 = np.array(df1)
            ar2 = []
            for i in ar1:
                x = i[0]
                y = i[1]
                l = y-x
                ix = i[3]
                iy = i[4]
                ibmean = (ix+iy)/2
                
                
                if l >= 5:
                    #For each insulation minima inside the large TAD ( in (x+1,y-1) ):
                    # -0.25 to the total scores of this TAD.
                    try:
                        minnumb = insmins[(insmins[:,0]>(x+1)) & (insmins[:,0]<(y-1))][:,1]
                        minnumb = minnumb[minnumb[:] < ibmean]
                        minnumb = -len(minnumb)*0.25
                    except IndexError:
                        minnumb = 0
                    except ValueError:
                        minnumb = 0
                else:
                    minnumb = 0
                #IMM = minimal insulation minimum inside the TAD.
                try:
                    imm = min(ff[(ff[:,0]>(x+1)) & (ff[:,0]<(y-1))][:,1])
                except IndexError:
                    imm = 26
                except ValueError:
                    imm = 26
                
                #If IMM is lower than the insulation value at each of the boundaries, 
                #then -3 to the total scores. 
                #...at one of the boundaries, then -1 to the total scores.
                if imm == 26:
                    inside = 0
                else:
                    if l >= 5:
                        if ((imm - ix)<0) & ((imm - iy)<0):
                            inside = -3
                        elif ((imm - ix)<0) & ((imm - iy)>=0):
                            inside = -1
                        elif ((imm - ix)>=0) & ((imm - iy)<0):
                            inside = -1
                        else:
                            inside = 0
                    else:
                        inside = 0

                try:
                    dlt_inside = max(deltas[(deltas[:,0]>(x+1)) & (deltas[:,0]<(y-1))][:,1])
                except IndexError:
                    dlt_inside = 0
                except ValueError:
                    dlt_inside = 0
                
                #For each strong delta minima (boundary strength >0.2)
                #inside the large TAD ( in (x+1,y-1) ):
                # -0.25 to the total scores of this TAD.
                if dlt_inside > 0.2:
                    d_insd = -9
                elif (dlt_inside >0) & (dlt_inside<=0.2):
                    d_insd = 0
                else:
                    d_insd = 0
                i = np.append(i,[inside,d_insd,minnumb]).tolist()
                ar2.append(i)
                
            arr = np.array(ar2)
            arr = np.hstack((arr,np.zeros((len(arr),1))))
            #See modul_score function description above.
            moduls = arr[:,2]
            p10 = np.percentile(moduls, 10)
            p25 = np.percentile(moduls, 25)
            p50 = np.percentile(moduls, 50)
            p75 = np.percentile(moduls, 75)
            p85 = np.percentile(moduls, 85)
            p90 = np.percentile(moduls, 90)
            p95 = np.percentile(moduls, 95)
            
            table = [
                [0, p10, -5],
                [p10, p25, -1],
                [p25, p50, 0],
                [p50, p75, 1],
                [p75, p85, 2],
                [p85, p90, 2.5],
                [p90, p95, 3],
                [p95, 1, 5],
            ]
            #SA = segment_aggregation parameter (modularity)
            arr[:,10] = modul_score(arr[:,2],table)
            dfg = pd.DataFrame(arr)
            updsegs.append(dfg)
        df3 = pd.concat(updsegs,axis=0,ignore_index = True)
        
        #Balancing large TADs.
        df3[10].loc[(df3[1]-df3[0])>8] *= 1.5
        df3[10].loc[(df3[1]-df3[0])>14] *= 1.25
        
        df3[11] = df3[[5,6,7,8,9,10]].sum(axis=1)
        #Merge with domains called from insulation minima.
        
        df3 = pd.concat([df3,dmins],axis=0,ignore_index = True)
        # Values of maximal insulation inside the domain
        # included in the column "SA_score"
        # for the domains called from insulation minima.
        
        df3[11].iloc[df3.loc[df3[11] == 0].index]=0.001
        df3[11].loc[(df3[1]-df3[0])<4] /= 2
        
        df3.columns = ['X','Y','Modul','Insul_x','Insul_y','IM_border','DLT_border',
                       'IMM_inside','DLT_inside','IM_inside','SA_score','Total']
#         print(df3)
        df3 = df3.sort_values('Total',ascending=False).reset_index(drop=True)
        df3 = df3.loc[df3.duplicated(subset=['X','Y'],keep='first') != True]

        df6 = df3.sort_values(['X','Y'], ascending=True).reset_index(drop=True)
        df6[['X','Y']] *= binsize
        df6 = np.array(df6)
        df6 = np.round(df6,5)
#         for i in df6:
#             i = i.tolist()
#             print(i)
        df3 = df3.sort_values(['Y'], ascending=True).reset_index(drop=True)
        df3 = df3[['X','Y','Total']]
        df3.columns = ['start','end','score']
        xylist = df3[['start','end']].values.tolist()
        xylist = [tuple(k) for k in xylist]
        scorelist =  df3['score'].values.tolist()
        listofcons,max_sc = consensus_domains(xylist,scorelist)
        
        consensus_df = pd.DataFrame(listofcons).sort_values([0,1], ascending=True).reset_index(drop=True)
        tad_len = np.mean(consensus_df[1] - consensus_df[0])
        
        consensus_df[0] = consensus_df[0]*binsize
        consensus_df[1] = consensus_df[1]*binsize-1
        consensus_df.insert(0, 0, chrom, allow_duplicates=True)
        consensus_df.columns = range(consensus_df.shape[1])
        consensus_df = consensus_df.astype({0:str,1:int,2:int})
        consensus_df = consensus_df[[0,1,2]]
        ann = consensus_df[[0,1,2,0,1,2]]
        mean_sc = max_sc/len(consensus_df)
        meta.append([sp,chrom,consensus_df,tad_len,mean_sc,g2])
        print(sp,chrom,g2,tad_len,mean_sc)
        if tad_len < sizes[-1]:
            print('bye')
            break
    print('end of chrom',sp,chrom)
    #     print ("{:g} s".format(time.process_time() - start_time))
meta = np.array(meta)

for sp, s in product(species,sizes):
    meta1 = meta[(meta[:,0]==sp)]

    meta4 = []
    for chrom in chroms:
        meta2 = meta1[(meta1[:,1]==chrom)]
        meta3 = meta2[np.argmin(abs(meta2[:,3]-s))]
        meta3 = meta3.tolist()
        meta4.append(meta3)
    meta4 = np.array(meta4)
    meta5 = meta4[:,2]
    bed = np.concatenate(meta5,axis=0)
    ann = np.concatenate([bed,bed],axis=1)
    
    path = '/home/al/Vp/0402/'
    try:
        makedirs(path)
    except FileExistsError:
        pass
    
    stat = meta4[:,[0,1,3,4,5]]
    stat1 = np.concatenate(stat,axis=0)
#     np.savetxt('/home/al/Vp/0402/'+sp+'_'+str(s)+'stat.bed',stat1, delimiter='\t', newline='\n',fmt='%s')
#     np.savetxt('/home/al/Vp/0402/'+sp+'_'+str(s)+'.bed',bed, delimiter='\t', newline='\n',fmt='%s')
#     np.savetxt('/home/al/Vp/0402/'+sp+'_'+str(s)+'.ann',ann, delimiter='\t', newline='\n',fmt='%s')
print ("{:g} s".format(time.process_time() - start_time))


# In[ ]:




