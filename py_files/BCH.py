#!/usr/bin/env python
# coding: utf-8

# ## BCH Code

# In[ ]:


get_ipython().run_line_magic('run', 'Conversions.ipynb')
get_ipython().run_line_magic('run', 'RS.ipynb')
#from sage.all import *
#from Conversions import _DetermineInput, _BitStringToInt, _IntToBitString, _IntToPol, _PolToInt
#from RS import *

class BCHCode:
    
    def __init__(self, n, b, D, q, shortening = 0):
        
        if q != 2:
            raise ValueError('works only in binary case')
            
        self.m = Mod(q,n).multiplicative_order()
        self.q = q
        self.b = b
        self.D = D
        self.d = D
        self.n = n
        
        if shortening.parent() == ZZ:
            self.shortening = shortening
        else:
            raise ValueError('The shortening input must be an integer')
        
        
        if (self.n != self.q**self.m - 1):
            raise ValueError('Invalid input values: n != q^m - 1')
        
        self.C_RS = RSCode(self.n, self.n - self.D + 1, self.q**self.m)
        
        # Initialize field
        self.F = GF(self.q) # base field
        self.EF = GF(self.q**self.m) # extension field
        self.R = PolynomialRing(self.F, 'x')
        self.x = self.R.gen()
        self.alpha = self.EF.primitive_element()
        
        # Constructing generator matrix
        self.cosets = self.cyclotomic_cosets(self.n, self.q, self.b, self.D)
        
        self.generator_poly = self.BCH_generator_polynomial(self.x, self.alpha, self.D, self.cosets)
        
        if not (self.generator_poly.divides(self.x**(self.n) - 1)):
            raise ValueError('generator_poly is not a generator polynomial')
            
        
        self.k = self.n - self.generator_poly.degree()
        
        self.G = matrix(self.F, self.k, self.n, lambda i,j : self.generator_poly[(j+((self.n)-i)) % (self.n)])
        self.G = self.G.echelon_form()
        self.n = self.n - self.shortening
        self.k = self.k - self.shortening
        
        
    def cyclotomic_cosets(self, n, q, b, D):
        # compute cyclotomic cosets
    
        cosets = []

        for i in range(b,b+D-1):
            coset = [(i * q**j) % n for j in range(0,n-1)]
            coset = list(set(coset))
            coset.sort()
            cosets.append(coset)
        return cosets

    def minimal_polynomial(self, coset, x, alpha):
        # compute minimal polynomial from one coset
        poly = 1
        for j in range(len(coset)):
            poly *= (x - alpha**coset[j])
        return poly

    def BCH_generator_polynomial(self, x, alpha, D, cosets):
        # compute generator polynomial
        poly = self.minimal_polynomial(cosets[0], x, alpha)
        for i in range(1,D-1):
            poly = LCM(poly,self.minimal_polynomial(cosets[i],x,alpha))
        
        return poly
    
    
    def Encoding(self, message, zeropad = True, out = 'bin', product_k = None):
        
        if not product_k:
            product_k = self.k + self.shortening
        else:
            product_k = (product_k // self.k) * (self.k + self.shortening)
        
        data_type = _DetermineInput(message, self.q)
        
        if data_type == 'pol':
            message = list(message)
        elif data_type == 'bin':
            message = list(message)
        elif data_type == 'int':
            pass
        else:
            raise ValueError('Wrong data type')
            
        
        rem = len(message) % self.k
        
        if rem != 0:
            if zeropad:
                message.extend([self.F(0)]*(self.k-rem))
            else:
                raise ValueError('k does not divide input size')      
        
        message_padded = []
        for i in range(0, len(message), self.k):
            message_padded.extend([self.F(0)]*self.shortening + message[i:i+self.k])

        c = []
        
        # Encoding each chunk of size k
        for i in range(0, len(message_padded), self.k + self.shortening):
            encoded_chunk = self.EncodeChunk(message_padded[i:i+self.k+self.shortening]) 
            c.extend(encoded_chunk[self.shortening:])
        
        c = vector(self.F, c)
        
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
        if len(chunk) != (self.k + self.shortening):
            raise ValueError('Invalid chunk size')
            
        c = vector(self.F, chunk) * self.G
        
        return c
                
    
    def Decoding(self, received, out = 'bin', product_n = None):
        
        if not product_n:
            product_n = self.n + self.shortening
        else:
            product_n = (product_n // self.n) * (self.n + self.shortening)
        
        data_type = _DetermineInput(received, self.q)
        
        if data_type == 'pol':
            pass
        elif data_type == 'int':
            received = _IntToPol(received, self.q)
        elif data_type == 'bin':
            received = _BitStringToInt(received, self.q)
            received = _IntToPol(received, self.q)
        else:
            raise ValueError('Wrong data type')
        
        
        received_padded = []
        for i in range(0, len(received), self.n):
            received_padded.extend([self.F(0)]*self.shortening + received[i:i+self.n].list())
        
        # Check input size
        if len(received_padded) % product_n != 0:
            raise ValueError('Invalid input size')
            
        d = []
        
        for i in range(0,len(received_padded),self.n+self.shortening):
            decoded_chunk = self.DecodeChunk(received_padded[i:i+self.n+self.shortening]) 
            d.extend(decoded_chunk[self.shortening:])
            
        d = vector(self.F, d)
            
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
            
    
    def DecodeChunk(self, chunk):
        
        if (len(chunk) != (self.n + self.shortening)):
            raise ValueError('Invalid input size')
        
        
        # Need to convert to extension field
        chunk = vector(self.EF, chunk)
          
        # Decode with RS decoder
        chunk = self.C_RS.DecodeChunk(chunk)
        
        chunk = self.C_RS.EncodeChunk(chunk)
        
        for i in range(len(chunk)):
            if (chunk[i] != self.F(0) and chunk[i] != self.F(1)):
                return vector(self.F, [0] * (self.k+self.shortening))
            
        chunk = chunk.Mod(self.q) 
        
        cols = self.G.pivots()
        G_independent = self.G.matrix_from_columns(cols)
        
        c = self.DecodeChunkBCH(chunk, cols, G_independent)
            
        return c
    
    def DecodeChunkBCH(self, chunk, cols, G_independent):
        
        chunk_independent = [chunk[i] for i in cols]
        
        return vector(self.F, chunk_independent) * G_independent.inverse()


# In[ ]:


#C = BCHCode(n = 2**10 - 1, b=1, D = 115, q = 2, shortening = 257) 
#C.n


# In[ ]:


#m = '00110011001100110110110111001000110101010101111'
#m = vector(GF(2), [1,1,0,0,1,1,0,0,1,1])
#m = [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1]
#m = ''
#for i in range(128):
#    e = ZZ.random_element(0,2)
#    m = m + str(e)
    
#print(m)

#m = '11001'
#c = C.Encoding(m, out = 'pol')
#c = C.Encoding(m, out = 'bin')
#c = C.Encoding(m, out = 'int')
#print('codeword: ', c)

#positions = []
#for i in range(31): # add delta errors
#    position = ZZ.random_element(0,C.n)
#    while position in positions:
#        position = ZZ.random_element(0,C.n)  
      
#    positions.append(position)
#    if c[position] == 1:
#        c[position] = 0
#    else:
#        c[position] = 1 # flip the bit

#d = C.Decoding(c)

#print('Decoding status: ', d == m)

#d = C.Decoding(c, out = 'pol')
#d = C.Decoding(c, out = 'int')
#d = C.Decoding(c, out = 'bin')
#print('decoded word: ', d)
#print('Decoding status: ', d == m)
#print(d)

