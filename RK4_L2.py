import matplotlib.pyplot as plt

# Read the file and extract the necessary information

var_values = []
with open('./output/partialconvergence.txt', 'r') as f:
    num_vars = int(f.readline().strip())
    time_steps = list(map(float, f.readline().strip().split()))
    
    for i in range(len(time_steps)):
        var_values.append(list(map(float, f.readline().strip().split())))

# Transpose the var_values array
var_values_transposed = list(zip(*var_values))

# Plot the variables along time
plt.plot(time_steps, var_values_transposed[0],  label='L2')

plt.xlabel('Time')
plt.ylabel(r'$\log((S_{9h} - S_{3h})/(S_{3h} - S_{h}))/\log(3)$')
plt.ylim(0,4.2)
plt.title('L2 Norm Variable Convergence Test')
plt.legend()
plt.savefig('python/Images/partial_conv.png', dpi=300)