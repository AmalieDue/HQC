#!/usr/bin/env python
# coding: utf-8

# ## Reed-Muller Code:

# In[ ]:


#get_ipython().run_line_magic('run', 'Conversions.ipynb')
from sage.all import *
from Conversions import _DetermineInput, _BitStringToInt, _IntToBitString, _IntToPol, _PolToInt

import numpy as np
import itertools

class RMCode:
    
    def __init__(self, r, m, q):
        
        if q != 2:
            raise ValueError('works only in binary case')
            
        self.r = r
        self.m = m
        self.q = q
        
        self.n = 2**self.m
        self.k = sum(binomial(self.m,i) for i in range(self.r + 1))
        self.d = 2**(self.m - self.r)
        
        # Initializing field
        self.F = FiniteField(self.q)
        
        self.G = matrix(self.F, self.k, self.n, self.generateG(self.r, self.m))
        
        
    def generateG(self, r, m):
        ones = np.repeat(1, 2**m)
        bases = [np.array([(i // 2**a) % 2 for i in range(2**m)]) for a in range(m)]
        if r == 1:
            G = np.concatenate([ones.reshape(1,2**m), np.stack(bases)])
            return G
        elif r > 1:
            others = []
            for i in range(2,r+1):
                others.extend([math.prod(x) for x in itertools.combinations(bases, i)])
            
            G = np.concatenate([ones.reshape(1,2**m), np.stack(bases), np.stack(others)])
            return G
        else:
            raise ValueError('wrong value of r')
    
    
    def Encoding(self, message, zeropad = True, out = 'bin', product_k = None):
        
        if not product_k:
            product_k = self.k
        else:
            product_k = product_k
        
        data_type = _DetermineInput(message, self.q)
        
        if data_type == 'pol':
            message = list(message)
        elif data_type == 'bin':
            message = list(message)
        elif data_type == 'int':
            pass
        else:
            raise ValueError('Wrong data type')
        
        rem = len(message) % product_k
        
        if rem != 0:
            if zeropad:
                message.extend([self.F(0)]*(product_k-rem))
            else:
                raise ValueError('k does not divide input size')    
                
        c = []
        
        # Encoding each chunk of size k
        for i in range(0, len(message), self.k):
            c.extend(self.EncodeChunk(message[i:i+self.k]))
            
        # convert c to 'pol'
        c = vector(self.F, c)
        
        # convert c if necessary
        if out == 'pol':
            return c
        elif out == 'int':
            c = _PolToInt(c, self.q)
            return c
        elif out == 'bin':
            c = _PolToInt(c, self.q)
            c = _IntToBitString(c, self.q)
            return c
        else:
            raise ValueError('Unrecognized output')
        
            
    def EncodeChunk(self, chunk): 
        # Encode a chunk of size k
        
        if len(chunk) != self.k:
            raise ValueError('Invalid chunk size')
            
        c = vector(self.F, chunk) * self.G

        return c
        
    
    def Decoding(self, received, out = 'bin', product_n = None):
        
        if not product_n:
            product_n = self.n
        else:
            product_n = product_n
        
        data_type = _DetermineInput(received, self.q)
        
        # convert received if necessary
        if data_type == 'pol':
            pass
        elif data_type == 'int':
            received = _IntToPol(received, self.q)
        elif data_type == 'bin':
            received = _BitStringToInt(received, self.q)
            received = _IntToPol(received, self.q)
        else:
            raise ValueError('Wrong data type')
        
        # Check input size
        if len(received) % product_n != 0:
            raise ValueError('Invalid input size')
            
        d = []
        
        for i in range(0,len(received),self.n):
            d.extend(self.DecodeChunk(received[i:i+self.n]))
            
        # convert d to pol
        d = vector(self.F, d)
            
        # convert d if necessary
        if out == 'pol':
            return d
        elif out == 'int':
            d = _PolToInt(d, self.q)
            return d
        elif out == 'bin':
            d = _PolToInt(d, self.q)
            d = _IntToBitString(d, self.q)
            return d
        else:
            raise ValueError('Unrecognized output')
                
            
    def DecodeChunk(self, word):
        # Reed Decoding algorithm
        
        if (len(word) != self.n):
            raise ValueError('Invalid chunk size')
        
        decoded = []
        
        variables = [i for i in range(1, self.m + 1)]
        r = self.r
        G = self.G
        
        while (r > 0):
            current_variables = []
            current_variables.append(list(itertools.combinations(variables,r)))
            current_variables = list(current_variables)
            
            for i in range(binomial(self.m,r)-1, -1, -1):

                # compute the complementary set T
                T = [i for i in range(1,self.m+1)]
                for j in range(r):
                    T.remove(current_variables[0][i][j])
                    
                combinations_T = self.DecodingConvert(T)
                combinations_S = self.DecodingConvert(list(current_variables[0][i]))
                votings = []
                for i in range(len(combinations_T)):
                    voting = 0
                    for j in range(len(combinations_S)):
                        voting += word[combinations_T[i] + combinations_S[j]]

                    votings.append(voting)
                decoded = self.MajorityDecoding(votings) + decoded
                
            G_part = matrix(self.F, binomial(self.m,r), self.n, G[-binomial(self.m,r):][:])
            word_part = vector(self.F, decoded[:binomial(self.m,r)]) * G_part
            word += word_part
            
            G = G[:-binomial(self.m,r)][:]
            
            r -= 1
                
        # decode the final (i.e. the first in the list because of reverse decoding) coefficient
        decoded = self.MajorityDecoding(word.list()) + decoded
        
        return vector(self.F, decoded)
    
    
    def DecodingConvert(self, elements): # takes as input one tuple
        combinations = []
        for i in range(len(elements)+1):
            combinations.append(list(itertools.combinations(elements,i)))
        combinations = list(combinations)
        integers = []
        for i in range(len(combinations)):
            for j in range(len(combinations[i])):
                l = [0 for i in range(self.m)]
                for k in range(len(combinations[i][j])):
                    l[combinations[i][j][k] - 1] = 1
                integers.append(ZZ(l,2))
        return integers

            
    def MajorityDecoding(self, word):
        # Decode by comparing the number of ones with the number of zeros
        
        if len(word) == 1:
            return word
        elif word.count(1) >= floor(len(word) / 2):
            return [1]
        elif word.count(1) < floor(len(word) / 2):
            return [0]
        else:
            return "decoding failure"


# In[ ]:


#C = RMCode(r =2, m = 4, q = 2)
#C.k


# In[ ]:


#m = vector(GF(2), [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0])
#m = [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1]
#m = '110011001'

#c = C.Encoding(m, out = 'pol')
#c = C.Encoding(m, out = 'int')
#c = C.Encoding(m, out = 'bin')
#print('codeword: ', c)

#d = C.Decoding(c, out = 'pol')
#d = C.Decoding(c, out = 'int')
#d = C.Decoding(c, out = 'bin')
#print('decoded word: ', d)


# In[ ]:





# In[ ]:




