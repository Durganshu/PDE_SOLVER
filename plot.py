import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
data = np.genfromtxt('reference_results.csv', delimiter=',')

data1 = data[1:,:]

L = 100
H = 100

row_values = range(0,L)

col = 0
itr = 0
data2 = np.zeros((L, H))

while(1):
    for row in row_values:
        data2[row,col] = data1[itr,2]
        itr = itr + 1
    col = col+1
    if(col == 100):
        break

    
plt.imshow(np.transpose(data2), extent=[-1, 1, -1, 1])
plt.colorbar()
plt.show()
plt.savefig("results.png")  #savefig, don't show
