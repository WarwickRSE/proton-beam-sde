import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta

# Read the data from fort.10
data = np.loadtxt('./Fortran/fort.20')

# Create histogram
fig, ax = plt.subplots(figsize=(10, 6))
counts, bins, patches = ax.hist(data, bins=50, edgecolor='black', density=True)

# Generate x values and beta PDF
x = np.linspace(0, 1, 100000)
y = beta.pdf(x, a=1, b=10000)

# Overplot beta function
ax.plot(x, y, 'r-', linewidth=2, label=f'Beta(α=1, β=5)')

ax.set_xlabel('Value')
ax.set_ylabel('Probability Density')
ax.set_title('Histogram of fort.20 with Beta Distribution')
ax.set_yscale('log')
ax.set_xlim(0, 0.01)
ax.set_ylim(1e-20, 10000)
ax.grid(True, alpha=0.3)
ax.legend()
plt.show()
