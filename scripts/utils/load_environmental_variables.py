#!/usr/bin/env python
# coding: utf-8

# In[1]:


from dotenv import load_dotenv, find_dotenv, dotenv_values
import os


# In[2]:


# find .env automatically by walking up directories until it's found
dotenv_path = find_dotenv()
# load up the entries as environment variables
load_dotenv(dotenv_path)

root_path = os.path.join(os.environ.get("ROOT_PATH"),'')
input_data_path = os.path.join(os.environ.get("INPUT_PATH"),'')
build_files_path = os.path.join(os.environ.get("BUILD_PATH"),'')
processed_data_path = os.path.join(os.environ.get("PROCESSED_PATH"),'')

n_cores = int(os.environ.get("N_CORES"))

