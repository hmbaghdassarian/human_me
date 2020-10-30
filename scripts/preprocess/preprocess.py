#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import zipfile


# In[ ]:


def unpack_files(data_files = '../data.zip', build_files = '../build_files.zip', 
                data_out = None, build_files_out = None):
    '''
    Unzips and moves input project files.
    
    data_file is the full/path/to/data.zip and build_files is full/path/to/build_files.zip. 
    Default is project root. 
    Optionally can provide output/path/directory_name for these. Default will unzip and place them in the same path.
    
    '''
    if not os.path.isfile(data_files):
        raise ValueError('data.zip does not exist at this path')
    if not os.path.isfile(build_files):
        raise ValueError('build_files.zip does not exist at this path')
    
    if data_out is None:
        data_out = data_files.split('.zip')[0]
    if build_files_out is None:
        build_files_out = build_files.split('.zip')[0]
    
    with zipfile.ZipFile(data_files, 'r') as zip_ref:
        zip_ref.extractall(data_out)
        
    with zipfile.ZipFile(build_files, 'r') as zip_ref:
        zip_ref.extractall(build_files_out)        
    
    
    return os.path.join(data_out), os.path.join(build_files_out)
    
    
def create_environment(input_data_path, build_files_path, root_path = '../',processed_data_path = '../',
                       n_cores = None):
    '''
    Creates a .env file to work with pydotenv. Input_data_path and build_files_path must be the output of preprocess.unpack_files.
   Root_path is the path/to/project_root/.
   Processed_data_path is the path/to/output_directory where output files will be stored. Default is project root.
    n_cores is an integer specifying the number of cores with which to parallelize certain processes.
    
    '''
    
    if not os.path.isdir(processed_data_path):
        os.mkdir(processed_data_path)
    
    if not os.path.isdir(input_data_path) or not  os.path.isdir(build_files_path):
        raise ValueError('Must specify data path and build_files path, see preprocess.unpack_files')
    
    processed_data_path, input_data_path, build_files_path, root_path = os.path.join(os.path.abspath(processed_data_path),''), os.path.join(os.path.abspath(input_data_path),''), os.path.join(os.path.abspath(build_files_path),''), os.path.join(os.path.abspath(root_path),'')
    
    if n_cores is not None:
        n_cores = str(round(n_cores))
    
    with open(root_path+'.env', 'w') as f:
        entry0 = 'ROOT_PATH="' + root_path + '"'
        entry1 = 'RAW_PATH="' + input_data_path + '"'
        entry2 = 'BUILD_FILES="' + build_files_path + '"'
        entry3 = 'PROCESSED_PATH="' + processed_data_path + '"'
        if n_cores is not None:
            entry4 = 'N_CORES="' + n_cores + '"'
        else:
            entry4 = ''
            
        f.write(entry0 + '\n' + entry1 + '\n' + entry2 + '\n' + entry3 + '\n' + entry4)    

