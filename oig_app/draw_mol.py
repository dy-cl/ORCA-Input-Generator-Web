import io
import base64
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw_mol(xyz_data):
    xyz_data = xyz_data.strip()

    lines = xyz_data.split('\n')

    symbols = []
    positions = []
    colors = []
    sizes = []

    color_dict = {'C': 'black', 'O': 'red', 'N':'blue', 'H': 'grey', 'S': 'yellow', 'F': 'lime', 'Cl': 'green', 'I': 'purple'}
    size_dict = {'C': 200, 'O': 300, 'N': 250, 'H': 50, 'S': 350, 'F': 100, 'Cl': 200, 'Br': 350, 'I': 400}

    bond_rules = {
        'C':  {'neighbours': 3, 'distance': 2.2},
        'O':  {'neighbours': 2, 'distance': 1.5},
        'N':  {'neighbours': 3, 'distance': 1.5},
        'H':  {'neighbours': 1, 'distance': 1.4},
        'S':  {'neighbours': 2, 'distance': 2.0},
        'F':  {'neighbours': 1, 'distance': 1.5},
        'Cl': {'neighbours': 1, 'distance': 2.0},
        'Br': {'neighbours': 1, 'distance': 2.1},
        'I':  {'neighbours': 1, 'distance': 2.2},
    }


    for line in lines:
        symbol, x, y, z = line.split()
        symbols.append(symbol)
        positions.append((float(x), float(y), float(z)))

        color = color_dict.get(symbol, 'purple')  #Default color is blue if symbol is not found in the dictionary
        size = size_dict.get(symbol, 100)  #Default size is 100 if symbol is not found in the dictionary

        colors.append(color)
        sizes.append(size)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for symbol, pos, color, size in zip(symbols, positions, colors, sizes):
        ax.scatter(pos[0], pos[1], pos[2], c=color, s=size)

    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            dist = sum((positions[i][k] - positions[j][k]) ** 2 for k in range(3)) ** 0.5
            if symbols[i] in bond_rules and symbols[j] in bond_rules:
                rule_i = bond_rules[symbols[i]]
                rule_j = bond_rules[symbols[j]]
                if dist < rule_i['distance'] and dist < rule_j['distance']:
                    ax.plot([positions[i][0], positions[j][0]],
                            [positions[i][1], positions[j][1]],
                            [positions[i][2], positions[j][2]], color='gray')
            else:
                if dist < 1.5:  #Default bonding condition for undefined atoms
                    ax.plot([positions[i][0], positions[j][0]],
                            [positions[i][1], positions[j][1]],
                            [positions[i][2], positions[j][2]], color='gray')

    ax.set_axis_off()

    image_stream = io.BytesIO()
    fig.savefig(image_stream, format='png')
    image_stream.seek(0)  
    
    return image_stream
   