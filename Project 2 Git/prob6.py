import matplotlib.pyplot as plt
import numpy as np

def fileextract(filename):
    v = [0]; v2 = [0]; v3 = [0]
    with open(filename,'r') as f:
        data = [line.strip() for line in f]
    N3 = int(len(data)/3) # 3 vectors of length len/3
    for i in range(N3):
        v.append(float(data[i]))
        v2.append(float(data[i+N3]))
        v3.append(float(data[i+2*N3]))
    v.append(0);v2.append(0);v3.append(0)
    return np.asarray(v), np.asarray(v2), np.asarray(v3)


v1, v2, v3 = fileextract('problem6a2.txt') # N = 10
v1b, v2b, v3b = fileextract('problem6b4.txt')# N = 100

n = len(v1); nb = len(v1b)
print('n: ',n)
h1 = 1/(n-1); h1b = 1/(nb-1)
h = 1/(n+1); hb = 1/(nb+1)
print(h*11)
a = -1/h**2; ab = -1/hb**2
d = 2/h**2; db = 2/hb**2
N = len(v1); Nb = len(v1b)
print('N: ',N)

xvec = np.zeros(len(v1)); xvecb = np.zeros(len(v1b))

lmbda = np.zeros(n)
v = np.zeros_like(lmbda)
x = np.zeros_like(lmbda)
vn = np.zeros_like(lmbda)

lmbdab = np.zeros(nb)
vb = np.zeros_like(lmbdab)
xb = np.zeros_like(lmbdab)
vnb = np.zeros_like(lmbdab)

for j in range(n):
    lmbda[j] = d + 2*a*np.cos((j*np.pi)/(N+1))
    v[j] = np.sin((j*np.pi)/(N+1))
    x[j] += j*h
for j in range(nb):
    lmbdab[j] = d + 2*a*np.cos((j*np.pi)/(Nb+1))
    vb[j] = np.sin((j*np.pi)/(Nb+1))
    xb[j] += j*hb
for i in range(n-1):
    xvec[i] += i*h1
for i in range(nb-1):
    xvecb[i] += i*h1b
#print(h1b*101)
xvec[-1] = 1; xvecb[-1] = 1
x = np.append(x,1) # Applying boundaries
v = np.append(v,0)
xb = np.append(xb,1); vb = np.append(vb,0)


fig, ax = plt.subplots(1,2,figsize=(10,6))
ax[0].set_title('Eigenvector numerical and analytic N = 10')
ax[1].set_title('N = 100')
ax[0].plot(xvec,v1/np.max(v1),label='lowest eig')
ax[0].plot(xvec,v2/np.max(v2),label='2nd lowest eig')
ax[0].plot(xvec,v3/np.max(v3),label='3rd lowest eig')
ax[0].plot(x,v,label='Analytical')
ax[1].plot(xvecb,v1b/np.max(v1b),label='lowest eig')
ax[1].plot(xvecb,v2b/np.max(v2b),label='2nd lowest eig')
ax[1].plot(xvecb,v3b/np.max(v3b),label='3rd lowest eig')
ax[1].plot(xb,vb,label='Analytical')
ax[0].set_xlabel('$\hat{x}_i$',fontsize=14); ax[0].set_ylabel('$v(\hat{x}_i)$',fontsize=14)
ax[1].set_xlabel('$\hat{x}_i$',fontsize=14)
ax[0].legend(loc=3); ax[1].legend(loc=8)
plt.tight_layout()
plt.show()
