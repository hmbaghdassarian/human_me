#!/usr/bin/env python
# coding: utf-8

import os

def abs_path(path):
    return os.path.join(os.path.abspath(path),'')

def create_environment(build_path, input_path = None, outdir = None,
                       n_cores = None):
    """Creates a .env file to work with pydotenv. 
    
    Parameters
    ----------
    build_path: str
        "full/path/to/build_dir" - see README to download appropriate files
    input_path: str, defaults to build path
        "full/path/to/inputs" - see README to download appropriate files
    outdir: str, defaults to project root
        "path/to/output_dir" - specifies where to store any output files (e.g., corrected M-model)
    n_cores: int
        number of cores to use for parallelization
    
    *Note: if downloading build and input files from human_me.preprocess.unpack_files function, 
    the input_path and build_path args must agree with the outputs of that function. 
    """
    root_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__), '../../'))
    
    if input_path is None:
        input_path = build_path
    if outdir is None:
        outdir = root_path
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    if not os.path.isdir(input_path) or not  os.path.isdir(build_path):
        raise ValueError('Must correctly specify path to input files and build files')
    
    
    outdir, input_path, build_path, root_path = abs_path(outdir), abs_path(input_path), abs_path(build_path),                                                 abs_path(root_path)
    
    if n_cores is not None:
        n_cores = str(round(n_cores))
    
    with open(root_path+'.env', 'w') as f:
        entry0 = 'ROOT_PATH="' + root_path + '"'
        entry1 = 'INPUT_PATH="' + input_path + '"'
        entry2 = 'BUILD_PATH="' + build_path + '"'
        entry3 = 'PROCESSED_PATH="' + outdir + '"'
        if n_cores is not None:
            entry4 = 'N_CORES="' + n_cores + '"'
        else:
            entry4 = ''
            
        f.write(entry0 + '\n' + entry1 + '\n' + entry2 + '\n' + entry3 + '\n' + entry4)    

