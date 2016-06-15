import os
import numpy as np
import nipy
import pandas as pd
import re
import gc
from config import *
import time
from sklearn.metrics.pairwise import euclidean_distances
from scipy import ndimage,stats
import sqlite3 as sql
#import seaborn as sns
from  matplotlib.pylab import plot, show, savefig, xlim, figure,hold, ylim, legend, boxplot, setp, axes,text, title

class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs



class Cluster(object):
    def __init__(self,name):
        self.name=name
        self.stat=None
        self.linked_samples=None
        self.expr=None

class Donor(object):

    def __init__(self,name):
        self.name=name
        self.probes=None
        self.expression=None
        self.anot=None
        self.pacall=None
        self.data_path=None
        self.data_type=None
        self.sample_map=None
        self.ontology=None

    def load(self):

        if os.path.isfile(os.path.join(DATA_DIR,self.name,'Probes.csv')):
            print 'Start reading {} csv data'.format(self.name)
            with Timer() as t:
                self.pacall=pd.read_csv(os.path.join(DATA_DIR,self.name,'PACall.csv'), header=None).as_matrix()[:,1:]
                self.probes=pd.read_csv(os.path.join(DATA_DIR,self.name,'Probes.csv'))
                self.expression=pd.read_csv(os.path.join(DATA_DIR,self.name,'MicroarrayExpression.csv'), header=None).as_matrix()[:,1:]
                self.anot=pd.read_csv(os.path.join(DATA_DIR,self.name,'SampleAnnot_edit.csv'))
            print 'Finished in {}s'.format(t.secs)
            self.data_type='csv'
            self.data_path=os.path.join(DATA_DIR,self.name)
            self.sample_map=nipy.load_image(SAMPLE_MRI[self.name])

        elif os.path.isfile(os.path.join(DATA_DIR,self.name,self.name+'.db')):
            print 'Connecting to {} sql database...'.format(self.name)
            con=sql.connect(os.path.join(DATA_DIR,self.name,self.name+'.db'))
            print 'Connected'
            self.sample_map=nipy.load_image(SAMPLE_MRI[self.name])
            self.data_path=os.path.join(DATA_DIR,self.name)
            self.data_type='sql'
        else:
            raise ValueError('There is no Allen Brain data in {}'.format(os.path.join(DATA_DIR,self.name)))


    def save_expression_map(self,gene_name, type='expr',samples_id=None):
        pass

    def get_gene_expr_info(self,gene_name,samples_id=None, probe_mode=None):

        if self.data_type=='csv':
            probe_ids=self.probes.gene_symbol.str.contains('^'+gene_name+'$')
            if np.sum(probe_ids.shape)==0:
                print 'Gene {} did not found!'.format(gene_name)
                return None

            probe_ids=np.where(probe_ids==1)[0]
            gene_expression=self.expression[probe_ids,:]
            gene_expression_bin=self.pacall[probe_ids,:]

            if probe_mode=='best' and np.sum(probe_ids)>1:
                cor = np.corrcoef(gene_expression)
                meancor = np.mean(cor, axis=0)
                maxrow = np.argmax(meancor, axis=0)
                gene_expression=gene_expression[maxrow:(maxrow+1),:]
                gene_expression_bin=gene_expression_bin[maxrow:(maxrow+1),:]

            elif probe_mode=='mean' and np.sum(probe_ids)>1: #TODO test
                gene_expression=gene_expression.mean(axis=1)
                gene_expression_bin=gene_expression_bin.mean(axis=1)

            elif probe_mode=='all':
                gene_expression= gene_expression
                gene_expression_bin= gene_expression_bin

            else:
                raise ValueError('Unknown probes mode')

            if samples_id is not None:
                mask=np.ones_like(range(gene_expression.shape[1]), dtype=bool)
                mask[samples_id]=False
                gene_expression_inside=gene_expression[:,mask]
                gene_expression_outside=gene_expression[:,~mask]
                gene_expression_inside_bin=gene_expression_bin[:,mask]
                gene_expression_outside_bin=gene_expression_bin[:,~mask]

                return {'inside':gene_expression_inside,  'outside':gene_expression_outside,
                       'inside_bin':gene_expression_inside_bin,  'outside_bin':gene_expression_outside_bin}

            return {'whole':gene_expression, 'whole_bin':gene_expression_bin}

        if self.data_type=='sql':
            raise ValueError('sql not implemented') #TODO

def form_clusters(data,threshold,type='p-value',cluster_size_threshold=1):

    s=ndimage.morphology.generate_binary_structure(3,3)
    if type=='p-value':
        clusters, n_clusters = ndimage.label((data < threshold) & (data>0),structure=s)
        stat_cl=ndimage.minimum(data,labels=clusters, index=range(1,n_clusters+1))
    elif type=='t-stat':
        clusters, n_clusters = ndimage.label(data > threshold,structure=s)
        stat_cl=ndimage.maximum(data,labels=clusters, index=range(1,n_clusters+1))
    else:
        raise ValueError('Wrong map type!')

    clusters_label=np.arange(1,n_clusters+1)
    count,sum=ndimage.measurements._stats(data,labels=clusters,index=clusters_label)
    clusters_mask=(count>cluster_size_threshold)
    if np.sum(count>10**5)!=0:
        raise ValueError('Some of the clusters are too huge for analysis {}.'
                         'Need to change the threshold to form clusters or check your input image.'
                         'If everything is correct, then probably you need to use -model correlation '.format(np.max(count))) #TODO correlation
    clusters_label=clusters_label[clusters_mask]
    return clusters,clusters_label,stat_cl


def link_samples2clusters(clusters_data,samples_data,dist_threshold=10):
    # return sample ids linked to clusters

    samples_mask=np.where(samples_data!=0)
    samples_N=len(samples_mask[0])
    sample_coordinate=np.array([np.array([samples_mask[0][i],samples_mask[1][i],samples_mask[2][i]]) for i in xrange(samples_N) ])

    N_clusters=clusters_data[1]
    linked_samples=[]
    cluster_statistic=[]

    for i in N_clusters:
        if i!=0:
            cluster_mask=np.where(clusters_data[0]==i)
            cluster_N=len(cluster_mask[0])
            cluster_coordinate=np.array([np.array([cluster_mask[0][j],cluster_mask[1][j],cluster_mask[2][j]]) for j in xrange(cluster_N) ])
            dist=euclidean_distances(sample_coordinate,cluster_coordinate)
            index=np.where(dist<=dist_threshold)
            index=np.unique(index[0])
            if index.shape[0]!=0:
                for l in index:
                    c=sample_coordinate[l,:]
                    s=samples_data[c[0],c[1],c[2]]
                    linked_samples.append(s)
                    #cluster_statistic.append(clusters_data[2][i-1])

    linked_samples=np.unique(np.array(linked_samples))

    return linked_samples if linked_samples.shape[0]!=0 else None



def sample_map_allen_space(image_path, annot_csv_path, save_path, type='well_id'):
    #assign values to samples in MNI space

    I=nipy.load_image(image_path)
    image_name=os.path.basename(image_path)
    df=pd.DataFrame.from_csv(annot_csv_path)
    coordinate, well_id=(np.array( df['mri_voxel_x']) , np.array(df['mri_voxel_y']), np.array(df['mri_voxel_z'] )), df[type]
    I._data[np.where(I._data!=0)]=0
    I._data[coordinate]=well_id
    nipy.save_image(I, os.path.join(save_path, image_name))


def plot_cluster_expression(out,data1,data2,donor,gene, image):

    # function for setting the colors of the box plots pairs
    def setBoxColors(bp):
        setp(bp['boxes'][0], color='blue')
        setp(bp['caps'][0], color='blue')
        setp(bp['caps'][1], color='blue')
        setp(bp['whiskers'][0], color='blue')
        setp(bp['whiskers'][1], color='blue')
        setp(bp['fliers'][0], color='blue')
        setp(bp['fliers'][1], color='blue')
        setp(bp['medians'][0], color='blue')

        setp(bp['boxes'][1], color='red')
        setp(bp['caps'][2], color='red')
        setp(bp['caps'][3], color='red')
        setp(bp['whiskers'][2], color='red')
        setp(bp['whiskers'][3], color='red')
        setp(bp['fliers'][2], color='red')
        setp(bp['fliers'][3], color='red')
        setp(bp['medians'][1], color='red')

    N_probes=data1.shape[0]

    fig = figure()
    ax = axes()
    hold(True)

    s=1
    f=2
    p_value=[]
    t_stat=[]
    ticks=[]
    for i in range(N_probes):
        t,p=stats.ttest_ind(data1[i,:], data2[i,:])
        bp = boxplot([data1[i,:],data2[i,:]], positions = [s, f], widths = 0.6)
        setBoxColors(bp)
        ticks.append( (s+f)/2. )
        s+=3
        f+=3
        p_value.append(p)
        t_stat.append(t)

    hB, = plot([1,1],'b-')
    hR, = plot([1,1],'r-')

    xlim(0,f+2)
    ylim(2,20)
    legend((hB, hR),('Inside', 'Outside'))

    for i in range(N_probes):
        text(f+3,10-i,'Probe #{}: p-value={}'.format(i+1,np.round(p_value[i],3) ) )
    ax.set_xticklabels(['probe #{}'.format(j) for j in range(1,N_probes+1)])
    ax.set_xticks(ticks)
    title('Donor {}, Allen Brain expression of gene {} inside/outside clusters formed in {} image'.format(donor, gene,image))

    hB.set_visible(False)
    hR.set_visible(False)
    try:
        savefig(os.path.join(out,donor+"_"+gene+".png") )
    except:
        savefig(os.path.join(out,donor+"_"+gene+".svg") )
