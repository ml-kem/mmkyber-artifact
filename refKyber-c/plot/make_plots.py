#   make_plots.py
#   === parse data and create some plots for the paper

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import FixedLocator, NullFormatter
from matplotlib.lines import Line2D
import sys

# Set font to Times New Roman
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.rm'] = 'Times New Roman'
#plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
#plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

args = ' '.join(sys.argv[1:])

if  args == 'speed':
    plot_speed = True
elif args == 'bytes':
    plot_speed = False
else:
    print("Pipe benchmark data in and invoke with argument 'speed' or 'bytes'. Thanks!")
    exit(0)

#   algorithm names, also plot labels
my_lab  = []
alg_lab = []
for sec in [ 128, 192, 256 ]:
    my_lab  +=  [   f'Kyber{4*sec}',
                    f'mmKyber-PKE-{sec}',
                    f'mmKyber-KEM-{sec}' ]
    alg_lab +=  [   f'ML-KEM-{4*sec}',  'mmKyber-PKE', 'mmKyber-KEM' ]

if plot_speed:
    my_var  = [ 'K-PKE.Encrypt()',
                'mmEnc()',
                'mmEncap()' ]
else:
    my_var  = [ 'ciphertext' ] * 3

# Data

x_values = [2**i for i in range(0, 11)]  # X-axis: 2^0 to 2^10

my_data  = [ [None] * len(x_values) for _ in my_lab ]

#   parse benchmark data from stdin

for line in sys.stdin:

    if line == '' or line[0] == '#':
        continue

    v = line.split()
    if len(v) < 6:
        continue

    if v[2] == 'N=':
        n = int(v[3])
        if n in x_values:
            x = x_values.index(n)
        else:
            continue
    else:
        continue

    if v[4] == 'len=':
        y = int(v[5])                       # bytes
    elif v[4] == 'cyc=':
        y = int(v[5]) / 1E6                 # millions of cycles
    else:
        continue

    for i in range(len(my_lab)):
        if v[0] == my_lab[i] and v[1] == my_var[i % 3]:
            my_data[i][x] = y



if plot_speed:
    #   make relative to kyber
    for i in range(0, 9, 3):
        for j in range(len(x_values)):
            my_data[i + 1][j] = my_data[i][j] / my_data[i + 1][j]
            my_data[i + 2][j] = my_data[i][j] / my_data[i + 2][j]
            my_data[i][j] = 1.0
else:
    #   divide by n
    for i in range(9):
        for j in range(len(x_values)):
            my_data[i][j] /= x_values[j]

#   baseline flat curve for reference
base1   = [ 1.0 for x in x_values ]

my_colour = [   'k', 'c', 'y',
                'k', 'g', 'm',
                'k', 'r', 'b' ]
my_marker = [   'o', 'D', 'v',
                'o', 'p', '<',
                'o', 's', '^' ]
# Plot
plt.figure(figsize=(12, 8))

for i in range(9):

    # don't take the baseline
    if plot_speed and i % 3 == 0:
        continue

    # mmKyber-KEM
    c = my_colour[i]
    m = my_marker[i]
    plt.plot(x_values, my_data[i], c + '-', linewidth=2, zorder=2)
    plt.plot(x_values, my_data[i], c + m, markersize=8, zorder=2)

# Kyber curve = 1.0

if plot_speed:
    plt.plot(x_values, base1, 'k-', linewidth=2, zorder=10)  # Increase z-order
    plt.plot(x_values, base1, 'ko', markersize=8, zorder=10)

    yoffs = [-4] * 9
    yoffs[0] = -8
    yoffs[4] = -10
    yoffs[5] = 2

    for i in [ 0, 1, 2, 4, 5, 7, 8 ]:
        plt.annotate(   f'${my_data[i][-1]:.2f} \\times$',
                        xy=(1, my_data[i][-1]), xytext=(8, yoffs[i]),
                        xycoords=('axes fraction', 'data'),
                        textcoords='offset points', color=my_colour[i],
                        fontsize=18)
else:
    for i in [ 0, 3, 6 ]:
        plt.annotate(   f'{alg_lab[i]}: ${my_data[i][-1]:.0f}$',
                        xy=(1, my_data[i][-1]), xytext=(8, -4),
                        xycoords=('axes fraction', 'data'),
                        textcoords='offset points', color=my_colour[i],
                        fontsize=18)

    for i in [ 7, 8 ]:
        plt.annotate(   f'{alg_lab[i]}: ${my_data[i-6][-1]:.0f}..{my_data[i][-1]:.0f}$',
                        xy=(1, my_data[i][-1]), xytext=(8, -4),
                        xycoords=('axes fraction', 'data'),
                        textcoords='offset points', color=my_colour[i],
                        fontsize=18)

# Set X-axis ticks
plt.xscale('log', base=2)

if not plot_speed:
    plt.yscale('log', base=2)

plt.gca().xaxis.set_major_locator(FixedLocator(x_values))  # Major ticks: 2^0, 2^2, ..., 2^10
plt.gca().xaxis.set_minor_locator(FixedLocator(x_values))  # Minor ticks: all 2^i
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"$2^{{{int(np.log2(x))}}}$"))  # Major tick labels
plt.gca().xaxis.set_minor_formatter(NullFormatter())  # Hide minor tick labels
plt.xticks(fontsize=20)  # Set font size

# Set Y-axis ticks
#plt.yticks([1, 10, 100, 1000], ["$10^0$", "$10^1$", "$10^2$", "$10^3$"], fontsize=16)
plt.yticks(fontsize=20)
# Disable minor ticks
#plt.tick_params(axis='y', which='minor', left=False)  # Disable minor ticks on Y-axis

# Custom legend
if plot_speed:
    legend_list = [
        Line2D( [0], [0], color='k', marker='o', linestyle='-',
                linewidth=3, markersize=10, label='Kyber/ML-KEM base') ]
else:
    legend_list = []

for i in range(0, 9):
    if i % 3 == 0:
        continue
    legend_list +=  [ Line2D(   [0], [0],
                        color=my_colour[i], marker=my_marker[i],
                        linestyle='-', linewidth=3, markersize=10,
                        label = my_lab[i])  ]

legend_list.reverse()

if plot_speed:
    plt.legend(handles=legend_list, loc='upper left', fontsize=16)
else:
    plt.legend(handles=legend_list, loc='lower left', fontsize=16)

# Add top and right borders
plt.gca().spines['top'].set_visible(True)
plt.gca().spines['top'].set_color('black')
plt.gca().spines['right'].set_visible(True)
plt.gca().spines['right'].set_color('black')

if plot_speed:
    plt.ylabel("Encryption cycles: Kyber / mmKyber", labelpad=10, fontsize=20)
else:
    plt.ylabel("Ciphertext bytes / Recipient", labelpad=10, fontsize=20)

plt.xlabel("Number of Recipients", labelpad=10, fontsize=20)

# Gridlines
plt.grid(False)

# Save as EPS file
if plot_speed:
    output_path = "speed_relative.pdf"
else:
    output_path = "bytes_relative.pdf"

plt.savefig(output_path, format='pdf', bbox_inches='tight')

# Show the plot
plt.tight_layout()

#plt.show()

print(f"Plot saved as {output_path}")
