#!/user/bin/env python3

"""
[Add module documentation here]

Author: Austin Egri
Date: 4/3/2019  

"""


import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde as kde



def submatsum(data,n,m):
    # return a matrix of shape (n,m)

    bs = data.shape[0]//n,data.shape[1]//m  # blocksize averaged over
    return np.reshape(np.array([np.sum(data[k1*bs[0]:(k1+1)*bs[0],k2*bs[1]:(k2+1)*bs[1]]) for k1 in range(n) for k2 in range(m)]),(n,m))

if __name__ == "__main__":
    # fname = sys.argv[1]

    # sample grid
    # size N
    
    '''
    # Uncomment to create sample text file
    N = 2**10 #2**15
    data = np.random.randint(0, high= 2, size= (N,N), dtype= bool)


    #print(data1.shape, data1)
    print(data.shape, data)
    np.savetxt("sample_data.txt", data, '%d', delimiter= ",")
    
    
    '''
    # load the data
    fname = "sample_data.txt"
    data = np.loadtxt(fname, delimiter= ",")
    #data = np.random.randint(0, high= 2, size= (2**15, 2**15), dtype= int)
    print(data)

    hist_data = submatsum(data, data.shape[0]//(2**5), data.shape[1]//(2**5))
    hist_data /= np.max(hist_data)


    #hist_data, edges = np.histogram(data)
    print(hist_data.shape)

    im = plt.imshow(hist_data, cmap=plt.cm.RdBu, interpolation='bilinear')
    plt.colorbar(im);

    plt.title("Heatmap of Game-of-Life Universe")

    #plt.hist(hist_data.shape)
    #sns.heatmap(kde(hist_data))

    plt.savefig("heatmap.jpg")
    #plt.show()





