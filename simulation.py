#!/usr/bin/env python
# coding: utf-8

# In[36]:

from sage.all import *
get_ipython().run_line_magic('run', 'HQC.ipynb')
import pickle
import time

w = 51
trials = 10000 // 8

current_time = time.strftime("%m.%d.%y %H:%M:%S", time.localtime())
output_name = 'out/concat RSRM1 w=%s %s.pkl' % (w, current_time)

# concat_RS_RM_1 = ConcatenatedCode(RSCode(n = 138, k = 16, q = 2**8), RMCode(r = 1, m = 7))


C = HQC(w = w , w_e = w , w_r = w, C = ConcatenatedCode(RSCode(n = 138, k = 16, q = 2**8), RMCode(r = 1, m = 7)), key_type = 'pol')
success = 0
minimum_decodable_error = C.n
maximum_decodable_error = 0
minimum_non_decodable_error = C.n
maximum_non_decodable_error = 0


for i in range(5):
    m = ''
    for i in range(C.k):
        e = ZZ.random_element(0,2)
        m = m + str(e)

    c = C.Encrypt(m, out = 'bin')
    if C.decodable_error_wt < minimum_decodable_error:
        minimum_decodable_error = C.decodable_error_wt
    if C.decodable_error_wt > maximum_decodable_error:
        maximum_decodable_error = C.decodable_error_wt
    if C.non_decodable_error_wt < minimum_non_decodable_error:
        minimum_non_decodable_error = C.non_decodable_error_wt
    if C.non_decodable_error_wt > maximum_non_decodable_error:
        maximum_non_decodable_error = C.non_decodable_error_wt

    d = C.Decrypt(c, out = 'bin')

    if d == m:
        success += 1

        
with open(output_name, "wb") as f:
    pickle.dump((C.num_ct_decoding_success, success, minimum_decodable_error, maximum_decodable_error, minimum_non_decodable_error, maximum_non_decodable_error),f)


# In[37]:


with open(output_name, "rb") as f:
    print(pickle.load(f))


# In[ ]:




