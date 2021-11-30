#!/usr/bin/env python
# coding: utf-8

# In[4]:

from sage.all import *
import time
import pickle
import RS.py

current_time = time.strftime("%m.%d.%y %H:%M:%S", time.localtime())
output_name = 'test %s.pkl' % (current_time)

F = GF(5)

with open(output_name, "wb") as f:
    pickle.dump([F(i) for i in range(10)] ,f)

with open(output_name, "rb") as f:
    print(pickle.load(f))


# In[ ]:





# In[ ]:




