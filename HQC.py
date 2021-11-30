#!/usr/bin/env python
# coding: utf-8

# ## HQC

# In[1]:


from sage.all import *
from ConcatenatedCode import *
from Conversions import *
#get_ipython().run_line_magic('run', 'ProductCode.ipynb')
#get_ipython().run_line_magic('run', 'ConcatenatedCode.ipynb')
#get_ipython().run_line_magic('run', 'Conversions.ipynb')


# In[11]:


class HQC:
    
    def __init__(self, w, w_e, w_r, C, key_type = 'pol', simulate_k = False):
        
        self.w = w
        self.w_e = w_e
        self.w_r = w_r
        self.C = C
        self.q = 2
        self.num_ct_decoding_success = 0
        self.simulate_k = simulate_k
        
        # if neither input nor output to code C are binary, then self.n and self.k are multiplied with the
        # power of the extension field size
        if (self.C.C1.q == self.C.C2.q and self.C.C1.q != self.q):
            self.n = self.C.n * self.C.C1.power
            self.k = self.C.k * self.C.C1.power
        else:
            self.n = self.C.n
            self.k = self.C.k
                    
        # initialize field and quotient ring
        self.F = GF(self.q)
        self.R = PolynomialRing(self.F, 'x'); x = self.R.gen()
        self.S = QuotientRing(self.R, x**self.n - 1, 'a'); a = self.S.gen()
        
        # generate public key and private key. The key data type is 'key_type', default is 'pol'.
        self.pk, self.sk = self.KeyGen(out = key_type)
        
        
    def KeyGen(self, out):
        # Input:
        # - out: Data type of output. Possible values are (default is 'pol'):
        #        'pol': Keys are elements of the ring S
        #        'int': Keys are lists of integers (only 1's and 0's)
        #        'bin': Keys are lists of bits
        
        # h is generated as a random element with random degree. Should there be a lower bound on the degree,
        # as we want the weight of h to be "big enough"?
        h = self.R.random_element(degree = (0, self.n))
        h = self.S(h)
        
        # x and y are randomly generated with wt(x) = wt(y) = self.w
        x = [0] * self.n
        x = self.RandomGenWithWeight(x, self.w)
        x = self.S(x)
        
        y = [0] * self.n
        y = self.RandomGenWithWeight(y, self.w)
        y = self.S(y)
        
        # compute syndrome
        s = x + h * y
        
        if self.simulate_k == True:
            # left part of parity check matrix, i.e. identity matrix
            H_left = matrix.identity(self.n)

            # right part of parity check matrix, i.e. rotation matrix of random h
            H_right = matrix(self.F, self.n, self.n, lambda i,j : h[(i - j) % self.n])

            # full parity check matrix
            H = block_matrix(1, 2, [H_left, H_left])

            self.k_random = 2 * self.n - rank(H)
            
        
        # convert keys if necessary
        if out == 'pol':
            pk = (h, s)
            sk = (x, y)
        elif out == 'int':
            pk = (_PolToInt(h, self.q), _PolToInt(s, self.q))
            sk = (_PolToInt(x, self.q), _PolToInt(y, self.q))
        elif out == 'bin':
            pk = (_PolToInt(h, self.q), _PolToInt(s, self.q))
            sk = (_PolToInt(x, self.q), _PolToInt(y, self.q))
            pk = (_IntToBitString(pk[0], self.q), _IntToBitString(pk[1], self.q))
            sk = (_IntToBitString(sk[0], self.q), _IntToBitString(sk[1], self.q))
        else:
            raise ValueError('Unrecognized output data type')
        
        return pk, sk
    
    
    def Encrypt(self, message, zeropad = True, out = 'bin'):
        # Input:
        # - message: Message data type error handling is handled by code C.
        # - out: Data type of output, i.e. ciphertext. Possible values are (default is 'bin'):
        #        'pol': Keys are elements of the ring S
        #        'int': Keys are lists of integers (only 1's and 0's)
        #        'bin': Keys are lists of bits 
        
        key_data_type = _DetermineInput(self.pk[0], self.q)
        
        # in order to make correct multiplication, the key data type should be 'pol'
        if key_data_type == 'pol':
            pk = self.pk
            sk = self.sk
        elif key_data_type == 'int':
            pk = (self.S(self.pk[0]), self.S(self.pk[1]))
            sk = (self.S(self.sk[0]), self.S(self.sk[1]))
        elif key_data_type == 'bin':
            pk = (_BitStringToInt(self.pk[0], self.q), _BitStringToInt(self.pk[1], self.q))
            pk = (self.S(pk[0]), self.S(pk[1]))
            sk = (_BitStringToInt(self.sk[0], self.q), _BitStringToInt(self.sk[1], self.q))
            sk = (self.S(sk[0]), self.S(sk[1]))
        
        
        c = self.C.Encoding(message, out = 'bin')
        # convert c to list of integers (1's and 0's) since list representation is necessary
        c = _BitStringToInt(c, self.q)
        
        u = []
        v = []
        
        non_decodable_error = [] # s * r2 + e'
        decodable_error = [] # x * r2 - r1 * y + e'
        
        for i in range(0, len(c), self.n):
            # e is randomly generated with wt(e) = self.w_e
            e = [0] * self.n
            e = self.RandomGenWithWeight(e, self.w_e)
            e = self.S(e)

            # r1 and r2 are randomly generated with wt(r1) = wt(r2) = self.w_r
            r1 = [0] * self.n
            r1 = self.RandomGenWithWeight(r1, self.w_r)
            r1 = self.S(r1)

            r2 = [0] * self.n
            r2 = self.RandomGenWithWeight(r2, self.w_r)
            r2 = self.S(r2)
            
            decodable_error_chunk = (self.sk[0] * r2 - r1 * self.sk[1] + e).list()
            decodable_error.extend(decodable_error_chunk)
            
            non_decodable_error_chunk = (pk[1] * r2 + e).list()
            non_decodable_error.extend(non_decodable_error_chunk)
            
            u.extend(r1 + pk[0] * r2)
            v.extend(self.S(c[i:i+self.n]) + pk[1] * r2 + e)
            
        self.decodable_error_wt = decodable_error.count(1)
        self.non_decodable_error_wt = non_decodable_error.count(1)
            
        # testing: Ensure that ciphertext cannot be decoded, i.e. ensure that wt(error) is too high
        if _DetermineInput(message, self.q) != 'bin':
            raise ValueError('testing is not going to work')
        
        v_testing = _PolToInt(v, self.q)
        v_testing = _IntToBitString(v_testing, self.q)
        
        v_decoded = self.C.Decoding(v_testing, out = 'bin')
        
        if v_decoded == message: # hopefully v_decoded and message are NOT equal
            self.num_ct_decoding_success = 1
    
        # convert ciphertext if necessary
        if out == 'pol':
            u, v = (self.S(u), self.S(v))
        elif out == 'int':
            u, v = (_PolToInt(u, self.q), _PolToInt(v, self.q))
        elif out == 'bin':
            u, v = (_PolToInt(u, self.q), _PolToInt(v, self.q))
            u, v = (_IntToBitString(u, self.q), _IntToBitString(v, self.q))
        else:
            raise ValueError('Unrecognized data type')
            
            
        return (u, v)
            
    
    def Decrypt(self, c, out = 'bin'):
        # Input:
        # - c: Ciphertext.
        # - out: Data type of output, i.e. plaintext. Possible values are (default is 'bin'):
        #        'pol': Keys are elements of the ring S
        #        'int': Keys are lists of integers (only 1's and 0's)
        #        'bin': Keys are lists of bits 
        
        
        # in order to make correct multiplication, the key data type should be 'pol'
        key_data_type = _DetermineInput(self.sk[0], self.q)
        
        # convert if necessary
        if key_data_type == 'pol':
            sk = self.sk
        elif key_data_type == 'int':
            sk = (self.S(self.sk[0]), self.S(self.sk[1]))
        elif key_data_type == 'bin':
            sk = (_BitStringToInt(self.sk[0], self.q), _BitStringToInt(self.sk[1], self.q))
            sk = (self.S(sk[0]), self.S(sk[1]))
            
            
        c_data_type = _DetermineInput(c[0], self.q)
        
        if c_data_type == 'pol':
            c = (_PolToInt(c[0], self.q), _PolToInt(c[1], self.q))
        elif c_data_type == 'int':
            pass
        elif c_data_type == 'bin':
            c = (_BitStringToInt(c[0], self.q), _BitStringToInt(c[1], self.q))
        
        
        decoded = []
        
        for i in range(0, len(c[0]), self.n):
            d = self.S(c[1][i:i+self.n]) - self.S(c[0][i:i+self.n]) * sk[1]
            d = _PolToInt(d, self.q)
            d = _IntToBitString(d, self.q)
            d = self.C.Decoding(d, out = 'bin')
            d = _BitStringToInt(d, self.q)
            decoded.extend(d)
        
        
        if out == 'pol':
            decoded = self.S(decoded)
        if out == 'int':
            pass
        elif out == 'bin':
            decoded = _IntToBitString(decoded, self.q)
        else:
            raise ValueError('Unrecognized data type')
        
        return decoded
        
        
    def RandomGenWithWeight(self, array, weight):

        for i in range(weight):
            cur = ZZ.random_element(0,len(array))
            while array[cur] == 1:
                cur = ZZ.random_element(0,len(array))
            
            array[cur] = 1
        
        return array


# In[ ]:


#C = HQC(w = 5, w_e = 5, w_r = 5, C = ProductCode(BCHCode(n = 15, b = 1, D = 7), RepetitionCode(n = 31)), key_type = 'bin')
#C = HQC(w = 7, w_e = 5, w_r = 5, C = ProductCode(RSCode(n=47, k=4, q=2**8), RSCode(n=47, k=4, q=2**8)), key_type = 'int')
#C = HQC(w = 10, w_e = 10, w_r = 10, C = ProductCode(RMCode(r = 1, m = 7), RMCode(r = 1, m = 7)), key_type = 'int')


#print('Public key: ', C.pk)
#print('\nSecret key: ', C.sk)
#C.k_random

