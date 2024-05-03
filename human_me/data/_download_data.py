import configparser
import sys
import os
import shutil

import gdown
from tqdm import tqdm

cur_dir = os.path.abspath(os.path.dirname(__file__))
def configure_data_path(data_path: str = cur_dir):
    """Configures the data path in human_me/data/data.ini.

    Parameters
    ----------
    data_path : str, optional
        path to store local large data files used for ME Model, by default os.path.abspath(os.path.dirname(__file__)
    """
    config = configparser.ConfigParser()
    config['data'] = {
        "local_data_path": data_path # defaults to human_me/data/
    }
    # TODO: add url data to config file instead of having a separate file_paths script?
    with open(os.path.join(cur_dir, 'data.ini'), 'w') as configfile:
        config.write(configfile)

def load_local_data_path():
    """Use the data.ini file to load the local data directory as a string"""
    if not os.path.isfile(os.path.join(cur_dir, 'data.ini')):
        raise NameError('The local data directory must be configured using "make build-data"')
    config = configparser.ConfigParser()
    config.read(os.path.join(cur_dir, 'data.ini'))

    local_dir = config['data']['local_data_path']
    build_local_path = os.path.join(local_dir, 'build/')
    input_local_path = os.path.join(local_dir, 'inputs/')
    return local_dir, build_local_path, input_local_path


google_files = {'build': {'psim_gold.h5': '1lpzOAhTEr_-4LGOmM4F2ucLLwBjGqJil', 
                           'recon2_2_only_psim.csv': '11JBC2rCRFhPdzN-nHn-QSdlCMD45iAS8'},
                'inputs': {'recon2_2.xml': '1Y5xdr5VtlJtA5r6OqicV_UNeRNhsOgcG'}
                }
prebuild_files = {'prebuild.zip': '1GCGHGi6sqgsnojL1-9avreOSx3ZAJA3y'}

# # tests
# google_files = {'build': {'test_1.txt': '1Qj9rQPqZzSlGB_utgKeFFKCQmA1qU_kB', 
#                                 'test_2.txt': '1Zw4LTLIm0of95nraeowwsvzrsDJO0MxU'}, 
#                     'inputs': {'test_3.txt': '1LuZGNd4h8PQ0FaF4n07ZJ88CWKlscIh5'}, 
#                     }
# prebuild_files = {'prebuild.zip': '1ReQ56JoxEyRDMi_e40_PNhoC6CazAXVO'}
def download_data(prebuild: bool = False):
    """Download the  large data files for ME Model building, and Recon2_2.xml as the default input M-model.

    Data should be downloaded somewhere with ~11 gb of space. 

    Parameters
    ----------
    prebuild : bool, optional
        whether to download "prebuild" files (used for generating the build files) in addition to building files, by default False
    """
    # inputs : bool, optional
    #     whether to download input metabolic model (Recon2.2xml) files in addition to building files, by default True
    local_dir, build_local_path, input_local_path = load_local_data_path()
    dir_map = {'build': build_local_path, 'inputs': input_local_path}

    # Make the data directories
    if os.path.exists(local_dir):
        if len(os.listdir(local_dir)) > 0:
            cont = input('The specified local data directory is not empty, do you wish to continue? Enter 1 for yes or 0 for no: ')
            if not bool(int(cont)):
                sys.exit('Please re-specify the local data directory must be configured using "make build-data"')
    else:
        os.makedirs(local_dir)
    for d_dir in [build_local_path, input_local_path]:
        if not os.path.exists(d_dir):
            os.makedirs(d_dir)
    
    # Download to the data directories
    for d_dir, file_dict in google_files.items():
        print('Download {} files'.format(d_dir))
        for file_name, file_id in file_dict.items():
            gdown.download(
                f"https://drive.google.com/uc?export=download&confirm=pbef&id={file_id}",
                output=os.path.join(dir_map[d_dir], file_name)
            )
#             gdown.download(id=file_id, output=os.path.join(dir_map[d_dir], file_name))
            # gdd.download_file_from_google_drive(file_id=file_id,
            #                                     dest_path=os.path.join(dir_map[d_dir], file_name),
            #                                     overwrite = False,
            #                                     showsize=True)
    if prebuild:
        print('Download prebuild files')
        prebuild_file = os.path.join(local_dir, 'prebuild.zip')
        gdown.download(id=prebuild_files['prebuild.zip'], output=prebuild_file)
        shutil.unpack_archive(prebuild_file, extract_dir=os.path.join(local_dir,'prebuild'))
        os.remove(prebuild_file)

if __name__ == "__main__":
    if sys.argv[1:]:
        configure_data_path(data_path = sys.argv[1]) # user provided
    else:
        configure_data_path() # default to data_path = human_me/human_me/data/
    if sys.argv[2:]:
        download_data(prebuild = bool(int(sys.argv[2])))
    else:
        download_data() # default to prebuild = False
