#!/usr/bin/env python
# coding: utf-8

# In[1]:


from dotenv import load_dotenv, find_dotenv, dotenv_values
import os


# In[2]:


print('load environmental variables')
# find .env automatically by walking up directories until it's found
dotenv_path = find_dotenv()
try:
    # load up the entries as environment variables
    load_dotenv(dotenv_path)

    root_path = os.environ.get("ROOT_DIR")
    print('The project root is: ' + root_path)
    local_data_path = os.environ.get("LOCAL_DATA_PATH")
    general_analyses_path = os.environ.get("GENERAL_ANALYSES_DIR")

    folder_id = os.environ.get("FOLDER_ID")
#     zip_name = os.environ.get("ZIP_NAME")

    email = os.environ.get("EMAIL")
    
except:
    raise ValueError('Create a .env file in root directory with appropriate variables filled in (see .env_template)')

