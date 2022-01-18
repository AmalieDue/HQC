#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('run', 'Conversions.ipynb')
#from sage.all import *
#from Conversions import _DetermineInput, _BitStringToInt, _IntToBitString, _IntToPol, _PolToInt

class RepetitionCode:
    
    def __init__(self, n, q):
        
        if q != 2:
            raise ValueError('works only in binary case')
            
        self.n = n
        self.k = 1
        self.d = n
        self.q = q
        
        self.F = GF(self.q)
        
        self.G = matrix(self.F, self.k, self.n, lambda i,j : 1)
        
    def Encoding(self, message, zeropad = True, out = 'bin', product_k = None):
        # message can be of the following types:
        # - 'pol': message = vector(GF(2), 10100110)
        # - 'int': message = [1,0,1,0,0,1,1,0]
        # - 'bin': messsage = '10100110'
        
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
        
        for i in range(len(message)):
            c.extend(self.EncodeChunk([message[i]]))
        
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
        elif data_type == 'bin':
            received = _BitStringToInt(received, self.q)
            received = _IntToPol(received, self.q)
        elif data_type == 'int':
            received = _IntToPol(received, self.q)
        else:
            raise ValueError('Wrong data type')
        
        # check input size
        if len(received) % product_n != 0:
            raise ValueError('Invalid input size')
            
        d = []
        
        for i in range(0,len(received),self.n):
            d.extend(self.DecodeChunk(received[i:i+self.n]))
            
        # convert d to 'pol'
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
        # decode chunk of size n
        
        if (len(word) != self.n):
            raise ValueError('Invalid chunk size')
        
        decoded = []
        
        # majority decoding
        if len(word) == 1:
            return vector(self.F, word)
        elif word.hamming_weight() >= floor(len(word) / 2):
            return vector(self.F, [1])
        elif word.hamming_weight() < floor(len(word) / 2):
            return vector(self.F, [0])
        else:
            return "decoding failure"


# In[ ]:


#C = RepetitionCode(5, 2)


# In[ ]:


#m = vector(GF(2), [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1])
#m = [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1]
#m = '11001100110011001'

#c = C.Encoding(m, out = 'pol')
#c = C.Encoding(m, out = 'int')
#c = C.Encoding(m, out = 'bin')
#print('codeword: ', c)

#d = C.Decoding(c, out = 'pol')
#d = C.Decoding(c, out = 'int')
#d = C.Decoding(c, out = 'bin')
#print('decoded word: ', d)


# In[ ]:




