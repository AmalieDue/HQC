#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#get_ipython().run_line_magic('run', 'HQC.ipynb')
import pickle
import time
from sage.all import *
from HQC import *

w = 77
trials = 500

current_time = time.strftime("%m.%d.%y %H:%M:%S", time.localtime())
output_name = 'output_files/SBCHRep1_w=%s/product SBCHRep1 w=%s %s.pkl' % (w, w, current_time)

# concat_RS_RM_1 = ConcatenatedCode(RSCode(n = 138, k = 16, q = 2**8), RMCode(r = 1, m = 7, q = 2))
# concat_RS_BCH_1 = ConcatenatedCode(RSCode(n = 139, k = 16, q = 2**8), BCHCode(n = 127, b = 1, D = 61, q = 2))
# product_SBCH_Rep_1 = ProductCode(BCHCode(n = 1023, b = 1, D = 115, q = 2, shortening = 257), RepetitionCode(n = 31, q = 2))

num_success = 0
minimum_decodable_error = 139 * 127
maximum_decodable_error = 0
minimum_non_decodable_error = 139 * 127
maximum_non_decodable_error = 0
num_ct_decoding_success = 0


for i in range(trials):
    C = HQC(w = w , w_e = w , w_r = w, C = ProductCode(BCHCode(n = 1023, b = 1, D = 115, q = 2, shortening = 257), RepetitionCode(n = 31, q = 2)), key_type = 'pol', simulation = True)
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
        
    if C.ct_decoding_success == 1:
        num_ct_decoding_success += 1

    d = C.Decrypt(c, out = 'bin')

    if d == m:
        num_success += 1

        
with open(output_name, "wb") as f:
    pickle.dump((num_ct_decoding_success, num_success, minimum_decodable_error, maximum_decodable_error, minimum_non_decodable_error, maximum_non_decodable_error),f)


# In[ ]:


#with open(output_name, "rb") as f:
#    print(pickle.load(f))

