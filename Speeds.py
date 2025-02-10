import matplotlib.pyplot as plt

# Read the data from the file
file_path = '/home/machado/Desktop/Universidade/Tese/thecode/output/Speeds.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()

# Process the data
first_line = lines[0].strip().split()
num_columns = len(first_line)
headers = [f'Column {i+1}' for i in range(num_columns)]
data = {header: [] for header in headers}

for line in lines:
    values = line.strip().split()
    for header, value in zip(headers, values):
        data[header].append(float(value))

# Plot the data
plt.figure(figsize=(10,8))
x = data[headers[0]]
plt.plot(x, data[headers[1]], label=r'$c_+$', color='black', linewidth=2)
plt.plot(x, data[headers[2]], label=r'$c_-$', color='red', linestyle='--', linewidth=2)
plt.plot(x, data[headers[3]], label=r'$c_+ + n_{cK}$', color='blue', linestyle='-.', linewidth=2)
plt.plot(x, data[headers[4]], label=r'$c_- + n_{cK}$', color='green', linestyle='-.', linewidth=2)
plt.plot(x, data[headers[5]], label=r'$c_+ + \lambda$', color='purple', linestyle=':', linewidth=3)
plt.plot(x, data[headers[6]], label=r'$c_- + \lambda$', color='orange', linewidth=2)
plt.plot(x, data[headers[8]], label=r'$c_+ + \mu$', color='brown', linestyle='-', linewidth=2)
plt.plot(x, data[headers[9]], label=r'$c_- + \mu$', color='green', linestyle='-', linewidth=2)
plt.xlim(0,1)
plt.axhline(y=0, color='gray', linewidth=1.5, label='0')
plt.xlabel('r')
plt.ylabel('Amplitude')
plt.title('Characteristic Speeds')
plt.legend()
plt.grid(False)
plt.tight_layout(pad=0.12)
#plt.savefig('python/Report_Thesis_Graphs/Speeds.pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.savefig('python/Report_Thesis_Graphs/Speeds_LG.pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)