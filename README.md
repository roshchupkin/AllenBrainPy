# AllenBrainPy
Python package to analyze Allen Humap Brain genes exression 

## Installation

Navigate to directory where you want to save this package and clone this repository:
     ```
     git clone https://github.com/roshchupkin/AllenBrainPy.git
     ```
## Update

You can update to the newest version using `git`. Navigate to your folder (where you cloned git repository):    
     ```
     git pull
     ```
     
## Usage

```
Python script to analyze VBM results and Allen Human Brain Atlas of gene
expression

optional arguments:
  -h, --help            show this help message and exit
  -o O                  path to save result folder
  -model {cluster_expression,correlation}
                        Analysis models
  -i I                  path input nifti image of VBM result map
  -d {all,caucasian}    choose all donors or only caucasian
                        (choices=['all','caucasian'])
  -threshold THRESHOLD  value threshold to form clusters from VBM result map
  -cl_size_threshold CL_SIZE_THRESHOLD
                        cluster size threshold to form clusters from VBM
                        result map
  -dist_threshold DIST_THRESHOLD
                        threshold for distance in voxels to link sample to
                        clusters
  -map_type {p-value,t-stat}
                        Type of VBM result map. For p-value clusters will be
                        formed for voxels < threshold,for t-stat clusters will
                        be formed for voxels > threshold, therefore first
                        split your image to negative and positive t-stat,save
                        as two maps with abs values and then run analysis
                        separately for both images.You can easily extrapolate
                        these two map types to any other result map.
  -result_name RESULT_NAME
                        name for saving results
  -probe_mode {all,best,mean}
                        gene name for expression analysis
  -plot                 plot boxplot of gene expression
  -gene_names GENE_NAMES [GENE_NAMES ...]
                        gene name for expression analysis
  -rsid RSID [RSID ...]
                        gene name for expression analysis

```