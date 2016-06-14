import os
import sys


DATA_DIR='/data/cluster_data/Allen_Brain/AB_data/'
MRI_DIR=''#TODO

DONOR=['donor10021',  'donor12876' , 'donor14380' , 'donor15496' , 'donor15697',  'donor9861']

DONOR_CAUCASIAN=['donor10021',  'donor12876', 'donor9861']

exclusions={'APOE4':'APOE', 'CELF1':'CUGBP1','FERMT2':'PLEKHC1', 'NME8':'TXNDC3'}

SAMPLE_MRI={
       'donor10021' :os.path.join(MRI_DIR,'rs_reorder_T1_donor10021.nii.gz'),
       'donor12876':os.path.join(MRI_DIR,'rs_reorder_T1_donor12876.nii.gz'),
        'donor14380':os.path.join(MRI_DIR,'rs_reorder_T1_donor14380.nii.gz'),
        'donor15496':os.path.join(MRI_DIR,'rs_reorder_T1_donor15496.nii.gz'),
        'donor15697':os.path.join(MRI_DIR,'rs_reorder_T1_donor15697.nii.gz'),
       'donor9861' :os.path.join(MRI_DIR,'rs_reorder_T1_donor9861.nii.gz')
            }


TEMPLATE="" #TODO define