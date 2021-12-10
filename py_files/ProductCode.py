#!/usr/bin/env python
# coding: utf-8

# ## Product Code

# In[83]:


#get_ipython().run_line_magic('run', 'RS.ipynb')
#get_ipython().run_line_magic('run', 'RM.ipynb')
#get_ipython().run_line_magic('run', 'BCH.ipynb')
#get_ipython().run_line_magic('run', 'Repetition.ipynb')
from sage.all import *
#from RS import *
#from RM import *
from BCH import *
from Repetition import *

import numpy as np

class ProductCode:
    
    def __init__ (self, code1, code2):
        
        self.C1 = code1
        self.C2 = code2
        
        self.n = self.C1.n * self.C2.n
        self.k = self.C1.k * self.C2.k
        self.d = self.C1.d * self.C2.d
        self.delta = floor( (self.C1.d - 1) / 2) * floor( (self.C2.d - 1) / 2)
        
        if self.C1.q != self.C2.q:
            raise ValueError('Wrong pair of codes: The fields are not the same')
            
    def Encoding(self, m, out = 'bin'):        
        
        c = self.C1.Encoding(m, out = 'int', product_k = self.k)
        
        c_total = []
        for i in range(0, len(c), self.C2.k * self.C1.n):
            c_chunk = np.reshape(c[i:i+self.C2.k * self.C1.n], (self.C2.k, self.C1.n))
            c_chunk = c_chunk.transpose()
            c_chunk = c_chunk.flatten()
            c_chunk = vector(ZZ, c_chunk)  
            c_total.extend(c_chunk)
                    
        c_total = self.C2.Encoding(c_total, out = out)
                
        return c_total
    
    
    def Decoding(self, r, out = 'bin'):
            
        d = self.C2.Decoding(r, out = 'int', product_n = self.n)
        
        d_total = []
        for i in range(0, len(d), self.C1.n * self.C2.k):
            d_chunk = np.reshape(d[i:i+self.C1.n * self.C2.k], (self.C1.n, self.C2.k))
            d_chunk = d_chunk.transpose()
            d_chunk = d_chunk.flatten()
            d_chunk = vector(ZZ, d_chunk)
            d_total.extend(d_chunk)
        
        d_total = self.C1.Decoding(d_total, out = out)
        
        return d_total


# In[87]:


#C = ProductCode(RSCode(n=47, k=4, q=2**8), RSCode(n=47, k=4, q=2**8))
#C = ProductCode(RSCode(n=127, k=8, q=2**8), RSCode(n=17, k=2, q=2**8))
#C = ProductCode(RMCode(r = 1, m = 3, q = 2), RMCode(r = 1, m = 3, q = 2))
#C = ProductCode(RMCode(r = 2, m = 8, q = 2), RMCode(r = 1, m = 6, q = 2))
#C = ProductCode(BCHCode(n = 1023, b = 1, D = 115, q = 2, shortening = 257), RepetitionCode(n = 31, q = 2))
#C.k


# In[89]:


#m = '101011'
#c = C.Encoding(m, out = 'bin')
#print("Codeword: ", len(c))

#m = '111011110101000111011110101000111011110101000111011110101000111011110101000111011110'

#print('Product code length: ', C.C1.n * C.C2.n)
#print('Product code delta: ', floor( (C.C1.d * C.C2.d - 1) / 2))
#m = [1,2,3,4,5,6,7]
#c = C.Encoding(m, out = 'bin')
#print(4*9*7)
#print("Codeword: ", len(c))
#d = C.Decoding(c, out = 'bin')
#print(d)

