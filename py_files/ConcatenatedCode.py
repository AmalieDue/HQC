#!/usr/bin/env python
# coding: utf-8

# In[6]:


#get_ipython().run_line_magic('run', 'RS.ipynb')
#get_ipython().run_line_magic('run', 'RM.ipynb')
#get_ipython().run_line_magic('run', 'BCH.ipynb')
from sage.all import *
from RS import *
from RM import * 

class ConcatenatedCode:
    
    def __init__ (self, code1, code2):
        
        self.C1 = code1
        self.C2 = code2
        
        self.n = self.C1.n * self.C2.n
        self.k = self.C1.k * self.C2.k
        self.d = self.C1.d * self.C2.d
        self.delta = floor( (self.C1.d - 1) / 2) * floor( (self.C2.d - 1) / 2)
        
        p, power = is_prime_power(self.C1.q, get_data = True)
        
        if self.C2.k != power:
            raise ValueError('Wrong pair of codes: Extension power of C1.q does not match C2.k')
        if self.C2.F != GF(p):
            raise ValueError('Wrong pair of codes: Extension field and base field do not fit together')
            
            
    def Encoding(self, m, out = 'bin'):
            
        c = self.C1.Encoding(m, out = 'bin')
    
        c = self.C2.Encoding(c, out = out)
        
        return c
    
        
    def Decoding(self, r, out = 'bin'):
        
        c = self.C2.Decoding(r, out = 'bin')
        
        c = self.C1.Decoding(c, out = out)
        
        return c


# In[7]:


#C = ConcatenatedCode(RSCode(n=138, k=16, q=2**8), RMCode(r = 1, m = 7, q = 2))
#C = ConcatenatedCode(RSCode(n=552, k=8, q=2**16), RMCode(r = 2, m = 5))
#C = ConcatenatedCode(RSCode(n=139, k=16, q=2**8), BCHCode(n = 127, b = 1, D = 61))


# In[8]:


#m = '10101010101011010'
#m = [1,2,3,4,5,6,0]
#c = C.Encoding(m, out = 'bin')
#print(c)
#d = C.Decoding(c, out = 'bin')
#print(d)


# In[ ]:




