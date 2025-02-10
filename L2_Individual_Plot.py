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
plt.plot(time_steps, var_values_transposed[0],  label='L2 Electric Field')
plt.plot(time_steps, var_values_transposed[1],  label='L2 Psi Field')
plt.plot(time_steps, var_values_transposed[2],  label='L2 Beta')
plt.plot(time_steps, var_values_transposed[3],  label='L2 Alpha')
plt.plot(time_steps, var_values_transposed[4],  label='L2 Phi')
plt.plot(time_steps, var_values_transposed[5],  label='L2 A')
plt.plot(time_steps, var_values_transposed[6],  label='L2 cphi')
plt.plot(time_steps, var_values_transposed[7],  label='L2 dphi')
plt.plot(time_steps, var_values_transposed[8],  label='L2 cPi')
plt.plot(time_steps, var_values_transposed[9],  label='L2 dPi')
plt.plot(time_steps, var_values_transposed[10], label='L2 trK')
plt.plot(time_steps, var_values_transposed[11], label='L2 Gammarr')
plt.plot(time_steps, var_values_transposed[12], label='L2 Chi')
plt.plot(time_steps, var_values_transposed[13], label='L2 Arr')
plt.plot(time_steps, var_values_transposed[14], label='L2 Lambdar')
plt.plot(time_steps, var_values_transposed[15], label='L2 Theta')

plt.xlabel('Time')
plt.ylabel(r'$\log((S_{9h} - S_{3h})/(S_{3h} - S_{h}))/\log(3)$')
#plt.ylim(3.8,4.2)
plt.title('L2 Norm Variable Convergence Test')
plt.legend()
plt.savefig('python/Images/partial_conv.png', dpi=300)