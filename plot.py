import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
data = np.genfromtxt('results.csv', delimiter=',')
plt.imshow(data, extent=[-1, 1, -1, 1])

plt.show()
plt.savefig("matplotlib.png")  #savefig, don't show
