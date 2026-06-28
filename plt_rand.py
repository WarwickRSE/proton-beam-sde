import numpy as np
import matplotlib.pyplot as plt

# Read the data from fort.10
data = np.loadtxt('./Fortran/fort.20')

# Create histogram
plt.figure(figsize=(10, 6))
plt.hist(data, bins=50, edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of fort.10')
plt.grid(True, alpha=0.3)
plt.show()
