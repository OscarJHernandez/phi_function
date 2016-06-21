import scipy.special 
import numpy as np
import math
from scipy import integrate

mpi = 139.9
hc =197.327


def G(x,q,r,nu,sigma,L):
	s0 = complex(0.0,1.0)*np.sqrt(0.25*q**2+mpi**2) 
	s1 = complex(0.0,1.0)*np.sqrt(0.25*(q**2)*(1.0-x**2)+mpi**2) 
	s2 = 0.5*q*x+s1
	#print s1,s2
	h = scipy.special.sph_jn(sigma,(s2*r/hc))[0][sigma]+complex(0.0,1.0)*scipy.special.sph_yn(sigma, (s2*r/hc))[0][sigma]
	dG = ((s2**(nu+1))/(s1))*h
	s = complex(0.0,1.0)*(np.pi/x)*dG*scipy.special.lpn(L, x)[0][L] 
	return np.real(s)

def phi(q,r,nu,sigma,L):
	s=integrate.quad(G, 0.0, 1.0,args=(q,r,nu,sigma,L,))[0]
	s = (1.0/(2.0*q))*s
	return s


nu =0
sigma = 0
L=0
q =0.01
x = 0.5
r=0.011
print x,q,G(x,q,r,nu,sigma,L)
print q,r,phi(q,r,nu,sigma,L),math.exp(-mpi*r/(hc))*(math.pi/(4.0*mpi))


nu =2
sigma = 0
L=0
q =0.1
x = 0.5
r=0.01
print x,q,G(x,q,r,nu,sigma,L)
print q,r,phi(q,r,nu,sigma,L)/(hc**2),math.exp(-mpi*r/(hc))*(math.pi/(4.0*r*hc))*(2.0-((mpi*r)/hc))


nu =2
sigma = 2
L=0
q =0.1
x = 0.5
r=0.01
print x,q,G(x,q,r,nu,sigma,L)
print q,r,phi(q,r,nu,sigma,L)/(hc**2),math.exp(-mpi*r/(hc))*(math.pi/(4.0*r*hc))*(1.0+((mpi*r)/hc))
