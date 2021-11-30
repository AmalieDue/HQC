#!/usr/bin/env python
# coding: utf-8

# ## Reed-Solomon code

# In[1]:


get_ipython().run_line_magic('run', 'Conversions.ipynb')

class RSCode:
    
    def __init__(self, n, k, q, alpha = None):
        
        if not (k < n and n <= q):
            raise ValueError('Invalid values for n, k, and q.')
            
        self.p0, self.power = is_prime_power(q, get_data = True)
            
        self.n = n
        self.k = k
        self.q = q
        
        self.d = self.n - self.k + 1
        
        self.tau = floor((self.n-self.k)/2)
        
        # Initializing field
        self.F = GF(self.q)
        self.R = PolynomialRing(self.F, 'X')
        self.p = self.F.primitive_element()
        
        # Constructing alpha-vector
        if not alpha:
            self.alpha = vector([self.p**i for i in range(self.n)])
        else:
            self.alpha = alpha
        
        # Constructing generator matrix
        self.G = matrix(self.F, k, n, lambda i,j : self.alpha[j]**i)
        
    def Encoding(self, m, zeropad = True, out = 'bin', product_k = None):
        
        if not product_k:
            product_k = self.k
        else:
            product_k = product_k
        
        # determine data type
        data_type = _DetermineInput(m, self.q)
        
        if data_type == 'int':
            m = _IntToPol(m, self.q)
        elif data_type == 'pol':
            pass
        elif data_type == 'bin':
            m = _BitStringToInt(m, self.q)
            m = _IntToPol(m, self.q)
        else:
            raise ValueError('Unrecognised input')
            
        
        rem = len(m) % product_k
        
        if rem != 0:
            if zeropad:
                m.extend([self.F(0)]*(product_k-rem))
            else:
                raise ValueError('k does not divide input size')
                
                
        c = []
        
        # Encoding each chunk of size k
        for i in range(0, len(m), self.k):
            c.extend(self.EncodeChunk(m[i:i+self.k]))
            
        c = vector(self.F, c)
        
        # Outputting decoded message in provided format
        if out == 'pol':
            return c
        elif out == 'bin':
            c = _PolToInt(c, self.q)
            return _IntToBitString(c, self.q)
        elif out == 'int':
            return(_PolToInt(c, self.q))
        else:
            raise ValueError('Error with output')
        
            
    def EncodeChunk(self, chunk):
        
        # Encode a chunk of size k
        if len(chunk) != self.k:
            raise ValueError('Invalid chunk size')
            
        c = vector(self.F, chunk) * self.G
        return c
    
    
    def Decoding(self, r, out = 'int', product_n = None):
        
        if not product_n:
            product_n = self.n
        else:
            product_n = product_n
        
        # Determining input data type
        data_type = _DetermineInput(r, self.q)

        if data_type == 'bin':
            r = _BitStringToInt(r, self.q)
            r = _IntToPol(r, self.q)
        elif data_type == 'int':
            r = _IntToPol(r, self.q)
        elif data_type == 'pol':
            pass
        else:
            raise ValueError('Unrecognised input')
        
        # Check input size
        if len(r) % product_n != 0:
            raise ValueError('Invalid input size')
            
        c = []
        
        for i in range(0,len(r),self.n):
            c.extend(self.DecodeChunk(r[i:i+self.n]))
            
        c = vector(self.F, c)
            
        
        # Outputting decoded message in provided format
        if out == 'pol':
            return c
        elif out == 'bin':
            tmp = _PolToInt(c, self.q)
            return _IntToBitString(tmp, self.q)
        elif out == 'int':
            return(_PolToInt(c, self.q))
        
        return c
    
    def DecodeChunk(self, chunk):
        
        if len(chunk) != self.n:
            raise ValueError('Invalid chunk size')
            
        # Constructing matrices
        M1 = matrix(self.F, self.n, self.tau + self.k, lambda i,j : self.alpha[i]**j)
        M2 = matrix(self.F, self.n, self.tau + 1, lambda i,j : chunk[i] * self.alpha[i]**j)
        M = M1.augment(M2)
        
        # Solving system
        RK = M.right_kernel()
        
        if len(RK.basis()) == 0:
            #return(None)
            return chunk[:self.k]
        
        sol = RK.basis()[0]

        # Constructing Q0 and Q1 polynomials
        Q0 = self.R(list(sol[:self.tau+self.k]))
        Q1 = self.R(list(sol[self.tau+self.k:]))

        # Calculating -Q0/Q1
        q, r = Q0.quo_rem(Q1)

        if r != 0:
            #print('Non-zero remainder (possibly >tau errors). Returning None')
            #return(None)
            return chunk[:self.k]

        out = []

        out.extend((-q).list())
        out.extend([self.F(0)]*(self.k-len(out)))

        return out


# In[4]:


#C = RSCode(9, 5, 2**4)


# In[1]:


#m = [1,0,1,0,1,0,1,1]
#m = [1,2,3,4,5,6,7]
#m = '0001001000110100010101100111'
#m = _IntToPol(m)
#m = [1, z4, z4 + 1, z4^2, z4^2 + 1, z4^2 + z4, z4^2 + z4 + 1, 0]
#print(m)


#c = C.Encoding(m, out = 'pol')
#c = C.Encoding(m, out = 'int')
#c = C.Encoding(m, out = 'bin')
#print('codeword: ', c)

#d = C.Decoding(c, out = 'pol')
#d = C.Decoding(c, out = 'int')
#d = C.Decoding(c, out = 'bin')
#print('decoded word: ', d)


# In[4]:


#m = []
#for i in range(C.k):
#    m.append(ZZ.random_element(0,C.q))
    
#print(m)
#delta = floor((C.d - 1)/2)
#print('delta = ', delta)
    
#c = C.Encoding(m)

#positions = []
#for i in range(delta): # add delta errors
#    position = ZZ.random_element(0,C.n)
#    while position in positions:
#        position = ZZ.random_element(0,C.n)  
      
#    positions.append(position)
#    c[position] = C.F(c[position] - 1) # flip the bit

#d = C.Decoding(c)

#print('Decoding status: ', d == m)


# In[ ]:




