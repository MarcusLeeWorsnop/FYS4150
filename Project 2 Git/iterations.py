import matplotlib.pyplot as plt
import numpy as np

iterations = np.array([2,10,7,31,36,64,81,112,138,183,199,247,270,336,366,432,485,557,600])
iterationsrand = np.array([2,8,14,30,42,63,87,115,141,170,207,236,289,338,388,436,514,549,629])
n = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])

plt.style.use('seaborn')
plt.title('Number of iterations of matrix N size',fontsize=14)
plt.plot(n,iterations,label='Tridiagonal')
plt.plot(n,iterationsrand,label='Dense matrix')
plt.xlabel('N matrix size',fontsize=14); plt.ylabel('Number of iterations',fontsize=14)
plt.legend(fontsize=12)
plt.yscale('log')
#plt.savefig('problem5b.pdf')
plt.show()
