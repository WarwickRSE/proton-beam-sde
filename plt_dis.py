import numpy as np
import matplotlib.pyplot as plt

# Read the data from fort.30
data = np.loadtxt('./Fortran/fort.20')

# Create x-axis from 0 to 1
x = np.linspace(0, 1, len(data))

# Create line plot
plt.figure(figsize=(10, 6))
plt.plot(x, data, linewidth=2)
plt.xlabel('x')
plt.ylabel('Value')
plt.title('Line plot of fort.30')
plt.xlim(0, 1)
plt.grid(True, alpha=0.3)
plt.show()
