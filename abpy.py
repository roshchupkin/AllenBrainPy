from tools import *
from config import *
import argparse
import gc



parser = argparse.ArgumentParser(description='Python script to analyze VBM results and Allen Human Brain Atlas of gene expression')

parser.add_argument("-o",required=True, type=str, help="path to save result folder")
parser.add_argument("-model",default='cluster_expression', choices=['cluster_expression','correlation'],type=str, help="Analysis models") #TODO extend
parser.add_argument("-i",required=True, type=str, help="path input nifti image of VBM result map")

group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument("-d",choices=['all','caucasian'], type=str, help="choose all donors or only caucasian (choices=['all','caucasian'])")
group1.add_argument("-donor",nargs='+', help="choose donor(s)")

parser.add_argument('-threshold',required=True,type=np.float64, help='value threshold to form clusters from VBM result map')
parser.add_argument('-cl_size_threshold',type=int,default=1, help='cluster size threshold to form clusters from VBM result map')
parser.add_argument('-dist_threshold',type=int,default=10, help='threshold for distance in voxels to link sample to clusters')
parser.add_argument('-map_type',choices=['p-value','t-stat'], required=True,
                    help='Type of VBM result map. '
                        'For p-value clusters will be formed for voxels < threshold,'
                        'for t-stat clusters will be formed for voxels > threshold, '
                        'therefore first split your image to negative and positive t-stat,'
                        'save as two maps with abs values and then run analysis separately for both images.'
                         'You can easily extrapolate these two map types to any other result map.')
parser.add_argument('-result_name', required=True, help='name for saving results')
parser.add_argument('-probe_mode', type=str,default='all',choices=['all','best', 'mean'], help='gene name for expression analysis')
parser.add_argument('-plot', action='store_true',default=False, help='plot boxplot of gene expression')

group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('-gene_names',nargs='+', help='gene name for expression analysis')
#group2.add_argument('-rsid',nargs='+', help='gene name for expression analysis') #TODO

args = parser.parse_args()
print args

I_map=nipy.load_image(args.i)

if args.d is not None:
    if args.d=='all':
        donors=DONOR
    elif args.d=='caucasian':
        donors=DONOR_CAUCASIAN
elif args.donor is not None:
    donors=args.donor
else:
    raise ValueError('one of the arguments -d -donor is required')

if args.model=='cluster_expression':

    cl=form_clusters(I_map._data,args.threshold,cluster_size_threshold=args.cl_size_threshold, type=args.map_type)

    with Timer() as t:
        I_map._data=cl[0]
        nipy.save_image(I_map,os.path.join(args.o,'cluster_label_'+str((args.threshold))+'_'+os.path.basename(args.i)))
    print 'Saved cluster label image {}'.format(os.path.join(args.o,'cluster_label_'+str((args.threshold))+'_'+os.path.basename(args.i)))

    gene_exression_info={}

    df_dic={}
    df_dic['p-value']=[]
    df_dic['t-stat']=[]
    df_dic['Donor']=[]
    df_dic['Gene']=[]
    df_dic['Linked Samples']=[]

    for d in donors:

        donor=Donor(d)
        donor.load()

        linked_samples=link_samples2clusters(cl,donor.sample_map._data,dist_threshold=args.dist_threshold)

        if not linked_samples is None:
            print 'Number of linked probes {}'.format(linked_samples.shape[0])
        else:
            print ("There is no linked samples to your map from donor {}".format(d))
            continue

        gene_exression_info[d]={}
        #gene_exression_info[d]['cluster_statistic']=cluster_statistic #TODO
        for g in args.gene_names:
            gene_exression_info[d][g]=donor.get_gene_expr_info(g,samples_id=linked_samples,probe_mode=args.probe_mode)
            N_probes=gene_exression_info[d][g]['inside'].shape[0]
            for i in range(N_probes):
                t,p=stats.ttest_ind(gene_exression_info[d][g]['inside'][i,:], gene_exression_info[d][g]['outside'][i,:])
                df_dic['p-value'].append(p)
                df_dic['t-stat'].append(t)
            df_dic['Donor']=df_dic['Donor']+[d]*N_probes
            df_dic['Gene']=df_dic['Gene']+[g]*N_probes
            df_dic['Linked Samples']=df_dic['Linked Samples']+[linked_samples.shape[0]]*N_probes
            if  args.plot:
                plot_cluster_expression(args.o,gene_exression_info[d][g]['inside'],gene_exression_info[d][g]['outside'],d,g,os.path.basename(args.i) )
        donor=None
        gc.collect()


    pd.DataFrame.from_dict(df_dic).to_csv(os.path.join(args.o,args.result_name +'.csv'))
    np.save(os.path.join(args.o,args.result_name+'.npy'),gene_exression_info)