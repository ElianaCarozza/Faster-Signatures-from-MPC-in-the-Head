import numpy as np
import math
from math import factorial
from math import gcd
from math import comb
from math import log

def phi(n):
    amount = 0        
    for k in range(1, n + 1):
        if gcd(n, k) == 1:
            amount += 1
    return amount

def L(blocksize):
	num = 0
	for k in range(1,blocksize+1):
		if k%2 == 1:
			D = gcd(blocksize-k,k)
			val = 0
			for d in range(1, D+1):
				if math.gcd(d,D)==d:
					val += phi(d)*factorial(int(blocksize/d))/(factorial(int((blocksize-k)/d))*factorial(int(k/d)))
			num += val
	return int(num/blocksize)

def part(n, k):
    def memoize(f):
        cache = [[[None] * n for j in range(k)] for i in range(n)]
        def wrapper(n, k, pre):
            if cache[n-1][k-1][pre-1] is None:
                cache[n-1][k-1][pre-1] = f(n, k, pre)
            return cache[n-1][k-1][pre-1]
        return wrapper
    @memoize
    def _part(n, k, pre):
        if n <= 0:
            return []
        if k == 1:
            if n <= pre:
                return [(n,)]
            return []
        ret = []
        for i in range(min(pre, n), 0, -1):
            ret += [(i,) + sub for sub in _part(n-i, k-1, i)]
        return ret
    return _part(n, k, n)

def max_numblocks(w,B):							# return the largest N such that u can have N distinct blocks yet generates less than B permutations
	N = 1
	numb = 1
	for j in range(w-N+2,w+1):
		numb *= j
	while numb <= B:
		N += 1
		numb = 1
		for j in range(w-N+2,w+1):
			numb *= j
	return N-1

def min_k(w,B):									# returns the smallest k (above w/2, by symmetry) such that c = (w choose k) is less than or equal to B
	k = int(w/2)
	c = math.comb(w,k)
	while c > B:
		k += 1
		c = math.comb(w,k)
	return (c,k)

def gvbound(w,bs,B,b=1,secpar=128):
	K = w*bs
	l = L(bs)
	N = max_numblocks(w,B)
	res = 0
	for n in range(2,N+1):						# starts at 2 to avoid the empty partition later when we isolate the biggest element
		val = 0
		(c,k) = min_k(w,B)						# computes the smallest k (above w/2, by symmetry) such that c = (w choose k) is less than or equal to B
		for i in range(k, w-n+2):				# for all possible values i of the first term in the partition of w (from k to w-(N-1))
			partitions = part(w-i,n-1)			# stores all partitions of the integer w-i into N-1 parts
			fact = 1
			for j in range(i+1, w+1):
				fact *= j						# computes w!/i!, a common factor of the terms summed over all partitions of w-i
			for partition in partitions:
				term = fact
				for v in partition:
					term = term//factorial(v) 	# computes w!/(i! * prod_v v!), where v goes over the elements of a partition of w-i
				if term <= B:					# checks that the term w!/(i! * prod_v v!) is below the bound B
					val += term					# val accumulates the terms of the sum accross all partitions of w-i, then accross all values of i
		res += factorial(n)*math.comb(l, n)*val 				# res accumulates the result accross all possible numbers n of blocks, from 1 to N = max_numblocks(w,B)
	dim = math.log((res+l)*(bs**w),2)			# dim is the base 2 logarithm of (K/w)^w * (the sum computed so far). /!\ will need to add lambda later
												# also adds l to the result because we excluded n=1
	return dim + b*secpar
