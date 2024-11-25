# plot_script.py
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_samples(data, output_path):
    # Convert the input data to a list if it is a tuple or array
    data = list(data)
    
    # Create a histogram of the data
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=30, color='blue', alpha=0.7)
    plt.title('Histogram of Generated Samples')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    
    # Save the plot to the provided output path
    plt.savefig(output_path)
    plt.close()

# This part is executed when called from the C++ code
def process_data_from_cpp(data_list, output_path):
    plot_samples(data_list, output_path)
