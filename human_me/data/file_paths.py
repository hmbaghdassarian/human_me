import os
from ._download_data import load_local_data_path

data_url = 'https://raw.githubusercontent.com/hmbaghdassarian/human_me_data/master/'
build_files_url = data_url + 'build/'
input_files_url = data_url + 'inputs/'

local_dir, build_local_path, input_local_path = load_local_data_path()
