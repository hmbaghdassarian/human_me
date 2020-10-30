#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pickle
import sys
import copy
import time

import cobra

import multiprocessing
import multiprocessing.pool
# from multiprocessing import Process
# from threading import Thread

sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
from utils import functions as func


# In[3]:


lp_path = '/data2/hratch/human_me/test_lp/'

with open(lp_path + 'me_model.pickle', 'rb') as handle:
    me_model = pickle.load(handle)

with open(lp_path + 'final_reactions.pickle', 'rb') as handle:
        final_reactions = pickle.load(handle)

print('Begin parallelization')


# In[14]:


def glpk(mu_val):
    print('Run glpk at growth {:.6f}'.format(mu_val))

    final_reactions_ = copy.deepcopy(final_reactions)
    for r in final_reactions_:
        if isinstance(r, func.ME_Reaction):
            if 'biomass' not in r.type:
                r.replace_coefficient_mu(mu_val = mu_val)
            else:
                r.replace_bound_mu(mu_val = mu_val, inplace = True)

    glpk_model = cobra.Model('glpk')
    glpk_model.add_reactions(final_reactions)
    glpk_model.objective = {glpk_model.reactions.biomass_dilution: 1}
    
    try:
        start = time.time()
        res = glpk_model.optimize()
        end = time.time()
        tot = str((end-start)/3600)
        
        with open(lp_path + 'optimization_outputs.tab', 'a') as f:
            f.write(str(mu_val) + '\t' + 'GLPK' + '\t' + tot + '\t' + '' + '\n')
        
        
        res = res.to_frame()
        res.to_csv(lp_path + 'glpk_' + str(mu_val).replace('.', '_') + '.csv')
    except:
        with open(lp_path + 'optimization_outputs.tab', 'a') as f:
            f.write(str(mu_val) + '\t' + 'GLPK' + '\t' + 'failed' + '\t' + 'failed' + '\n')
    
    del final_reactions_
    del glpk_model
    del res

def qminos(mu_val, precision):
    start = time.time()
    xq,statq,hsq = me_model.solve_lp(mu_val = mu_val, precision = precision)
    end = time.time()
    tot = str((end-start)/3600)
    
    with open(lp_path + 'optimization_outputs.tab', 'a') as f:
        f.write(str(mu_val) + '\t' + 'QMINOS_' + precision + '\t' + tot + '\t' + str(statq.max()) + '\n')
    
    res = pd.DataFrame(xq)
    res.to_csv(lp_path + 'qminos_' + precision + '_' + str(mu_val).replace('.', '_') + '.csv')
    del res
    del hsq

def qminos_double(mu_val):
    print('Run qminos double at growth {:.6f}'.format(mu_val))
    qminos(mu_val, precision = 'double')
    
def qminos_quad(mu_val):
    print('Run qminos quad at growth {:.6f}'.format(mu_val))
    qminos(mu_val, precision = 'quad')


# In[ ]:


#https://stackoverflow.com/questions/2957116/make-2-functions-run-at-the-same-time - threading
#https://stackoverflow.com/questions/7207309/how-to-run-functions-in-parallel - parallelizing
#https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic***

# def runSimultaneously(mu_val): #,*fns):
#     proc = []
#     fns = [glpk, qminos_double, qminos_quad]
#     for fn in fns:
#         p = Process(target=fn, args = (mu_val,))
#         p.start()
#         proc.append(p)
#     for p in proc:
#         p.join()   

# def runSimultaneously(mu_val):
#     Thread(target=glpk, args = (mu_val,)).start()
#     Thread(target=qminos_double, args = (mu_val,)).start()
#     Thread(target=qminos_quad, args = (mu_val,)).start()

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def run(mu_val, fctn):
    fctn(mu_val)
    
def runSimultaneously(mu_val):
    pool = multiprocessing.Pool(2)
    pool.starmap(run, zip([mu_val]*3, [qminos_double, qminos_quad]))#, glpk]))
    pool.close()
    pool.join()


# In[ ]:


with open(lp_path + 'optimization_outputs.tab', 'w') as f:
    f.write('mu' + '\t' + 'Algorithm' + '\t' + 'Run Time' + '\t' + 'Status' + '\n')

mu_vals = [1e-9, 1e-6, 1e-3, 5e-3, 1e-2, 2e-2, 3e-2, 5e-2, 1]
# mu_vals = # input a new list and wil optimization_outputs.tab will continue to be appended 
n_cores = len(mu_vals)

start = time.time()
pool = MyPool(3, maxtasksperchild=500) # run before to save memory
pool.map(runSimultaneously, mu_vals)        
pool.close()
pool.join()
end = time.time()
print('Total run time: {} (hrs)'.format((end-start)/3600))

