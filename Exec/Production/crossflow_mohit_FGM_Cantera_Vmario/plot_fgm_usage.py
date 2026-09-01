import yt
import matplotlib.pyplot as plt

# 1. Load the PeleLMeX plotfile (point this to your 'plt' folder)
ds = yt.load("plt00951") 
ad = ds.all_data()

# 2. Extract the variables that define your FGM manifold
# Based on your Header, we will use ZMIX and progress_variable
zmix = ad[("boxlib", "Y(ZMIX)")]
prog = ad[("boxlib", "progress_variable")]

# Extract a scalar to color the points (e.g., Temperature from the manifold)
temp = ad[("boxlib", "MANI_T")] 

# 3. Create the scatter plot
plt.figure(figsize=(8, 6))
# Using s=1 for small point sizes, alpha=0.5 for transparency if there are many overlapping points
scatter = plt.scatter(zmix, prog, c=temp, s=1, alpha=1.0, cmap='inferno')

# 4. Format the plot
plt.xlabel('Mixture Fraction [Y(ZMIX)]')
plt.ylabel('Progress Variable [progress_variable]')
#plt.title('Simulation Points in FGM Phase Space')
plt.colorbar(scatter, label='Temperature (K) [MANI_T]')

# Set limits if you know the bounds of your FGM table (e.g., 0 to 1)
# plt.xlim(0, 1.0)
# plt.ylim(0, 1.0)

plt.tight_layout()
plt.show()