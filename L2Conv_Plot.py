import matplotlib.pyplot as plt

# Read the file and extract the necessary information
with open('./output/convergence.txt', 'r') as f:
    num_vars = int(f.readline().strip())
    time_steps = list(map(float, f.readline().strip().split()))
    var_values = list(map(float, f.readline().strip().split()))

# Plot the variables along time
plt.plot(time_steps, var_values, label='log L2 Norm')

plt.xlabel('Time')
plt.ylabel(r'$\log((S_{4h} - S_{2h})/(S_{2h} - S_{h}))/log(4)$')
plt.ylim(0,5)
plt.title('L2 Norm Convergence Test')
plt.legend()
plt.savefig('python/Images/conv.png', dpi=300)