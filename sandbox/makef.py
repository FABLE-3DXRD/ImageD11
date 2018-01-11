
from __future__ import print_function



def one_item_func( grain, par ):
    diff = grains[i] - par*42. - 1.
    return diff*diff

def make_sum_func( n ):
    call = 'lambda p_'+', p_'.join([str(i) for i in range(n)]) 
    call += ' : '
    for i in range(n):
       call +=  '+ one_item_func( %d, p_%d )'%(i,i)
    return eval( call )

f = make_sum_func(3)
grains = {}
for i in range(10):
    grains[i] = 43

import minuit
m = minuit.Minuit( f )

m.migrad()
print(m.values)

f = make_sum_func(4)
grains = {}
for i in range(4):
    grains[i] = -3.

print('HELLO')
m = minuit.Minuit( f )

m.migrad()
print(m.values)