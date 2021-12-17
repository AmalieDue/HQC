#!/usr/bin/env python
# coding: utf-8

# ## Conversions

from sage.all import *

# In[19]:


def _DetermineInput(data, q):
    # determine data type
    F = GF(q)

    if isinstance(data, str) and all([(bit == '1' or bit == '0') for bit in data]):

        return('bin')

    if data[0] == None:
        return('none')
    elif data[0].parent() == F:
        return('pol')
    elif data[0].parent() == ZZ:
        return('int')
    else:
        return('unknown')
        

def _PolToInt(pol_array, q):

    # convert array of polynomial representation to array of integers
    
    F = GF(q)
    p, power = is_prime_power(q, get_data = True)

    pol_out = []

    for pol in pol_array:
        if not pol in F:
            raise ValueError('Invalid symbol')

        pol_out.append(ZZ(pol.polynomial().coefficients(sparse = False), base = p))


    return pol_out


def _IntToBitString(int_array, q):
    # Converts array of integers less than 2 to bit string
    
    p, power = is_prime_power(q, get_data = True)

    if any([(item > (q - 1) or item < 0) for item in int_array]):
        raise ValueError('Invalid integer values')

    number_of_bits = '0' + str(power) + 'b'

    return(''.join([format(item, number_of_bits) for item in int_array]))


def _BitStringToInt(bit_string, q):
    
    p, power = is_prime_power(q, get_data = True)
    F = GF(q)

    #if(self.C1.q != 2**self.C2.k):
    #    raise ValueError('Invalid field size for byte representation')

    # Converts array of 8 bit binary representations to integers

    #if (len(byte_string) % self.k_RM != 0) or any([not (bit == '1' or bit == '0') for bit in byte_string]):
    #    raise ValueError('Invalid byte string')

    #bit_string = list(bit_string)
    
    rem = len(bit_string) % power
    
    if rem != 0:
        bit_string += '0' * (power - rem)
        
    
    m_out = []
    for i in range(0,len(bit_string), power):
        current = bit_string[i:i+power]
        m_out.append(int("".join(str(x) for x in current), p))

    return m_out


def _IntToPol(array, q):
    # Convert array of integers less than q to elements of field

    F = GF(q)
    p, power = is_prime_power(q, get_data = True)
    
    if q == 2:
        return vector(F, array)
    elif q > 2:
        array_out = []

        for i in array:
            if not i < q:
                raise ValueError('Invalid symbol')
            array_out.append(F(ZZ(i).digits(p)))

        return array_out
    else:
        raise ValueError('Can not convert from list of integers to list of elements from field')
        
    
    


# In[ ]:




