#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#get_ipython().run_line_magic('run', 'HQC.ipynb')
import pickle
import time
from sage.all import *
from HQC import *

w = 62
trials = 1000

current_time = time.strftime("%m.%d.%y %H:%M:%S", time.localtime())
output_name = 'output_files/RMRep2_w=%s_new_h/product RMRep2 w=%s %s.pkl' % (w, w, current_time)

# concat_RS_RM_1 = ConcatenatedCode(RSCode(n = 138, k = 16, q = 2**8), RMCode(r = 1, m = 7, q = 2))
# concat_RS_BCH_1 = ConcatenatedCode(RSCode(n = 139, k = 16, q = 2**8), BCHCode(n = 127, b = 1, D = 61, q = 2))
# product_SBCH_Rep_1 = ProductCode(BCHCode(n = 1023, b = 1, D = 115, q = 2, shortening = 257), RepetitionCode(n = 31, q = 2))
# single_RM_1 = RMCode(r = 3, m = 14, q = 2)
# single_RM_2 = RMCode(r = 2, m = 15, q = 2)
# product_RM_RM_1 = ProductCode(RMCode(r = 2, m = 7, q = 2), RMCode(r = 1, m = 8, q = 2))
# product_RM_RM_2 = ProductCode(RMCode(r = 1, m = 6, q = 2), RMCode(r = 2, m = 8, q = 2))
# product_RS_RS_1 = ProductCode(RSCode(n = 47, k = 4, q = 2**8), RSCode(n = 47, k = 4, q = 2**8))
# product_SBCH_BCH_1 = ProductCode(BCHCode(n = 2**8 - 1, b = 1, D = 31, q = 2, shortening = 116), BCHCode(n = 127, b = 1, D = 57, q = 2))
# product_SBCH_RM_1 = ProductCode(BCHCode(n = 255, b = 1, D = 31, q = 2, shortening = 117), RMCode(r = 1, m = 7, q = 2))
# product_RM_Rep_1 = ProductCode(RMCode(r = 4, m = 9, q = 2), RepetitionCode(n = 31, q = 2))
# product_RM_Rep_2 = ProductCode(RMCode(r = 4, m = 9, q = 2), RepetitionCode(n = 33, q = 2))


#print(time.localtime())
num_success = 0
minimum_decodable_error = 2**9 * 33
maximum_decodable_error = 0
minimum_non_decodable_error = 2**9 * 33
maximum_non_decodable_error = 0
num_ct_decoding_success = 0
#low_h_wt = []
minimum_weight_h = 2**9 * 33
maximum_weight_h = 0



for i in range(trials):
    C = HQC(w = w, w_e = w, w_r = w, C = ProductCode(RMCode(r = 4, m = 9, q = 2), RepetitionCode(n = 33, q = 2)), key_type = 'pol', simulation = True, single_code = False)
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
    if C.w_h < minimum_weight_h:
        minimum_weight_h = C.w_h
    if C.w_h > maximum_weight_h:
        maximum_weight_h = C.w_h
    
        
    if C.ct_decoding_success == 1:
        num_ct_decoding_success += 1
        #low_h_wt.append(C.low_h_wt)

    d = C.Decrypt(c, out = 'bin')

    if d == m:
        num_success += 1

#print(time.localtime())
with open(output_name, "wb") as f:
    pickle.dump((trials, minimum_weight_h, maximum_weight_h, num_ct_decoding_success, num_success, minimum_decodable_error, maximum_decodable_error, minimum_non_decodable_error, maximum_non_decodable_error),f)


# In[ ]:


#print(time.localtime())
#with open(output_name, "rb") as f:
#    print(pickle.load(f))


# In[ ]:


#3, 14: 36,37
#2, 15: 83,84

