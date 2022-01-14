import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True

# Enter 'results.csv' for plotting numerical solution
# Enter 'unit_test_results.csv' for plotting analytical solution
filename = '../results/results.csv'
data = np.genfromtxt(filename, delimiter=',')

data1 = data[1:,:]

L = 101
H = 101

row_values = range(0,101)

col = 0
itr = 0
data2 = np.zeros((101, 101))


column = 2
if(filename == 'unit_test_results.csv'):
    column = 3

    
while(1):
    for row in row_values:
        data2[row,col] = data1[itr,column]
        itr = itr + 1
    col = col+1
    if(col == 101):
        break

    
plt.imshow(np.transpose(data2),cmap = cm.jet, extent=[0, 1, 0, 1])
plt.colorbar()
plt.show()

plt.savefig("results.png")  #savefig, don't show
