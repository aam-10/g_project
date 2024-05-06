import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
import os
import random
import statistics

colors = {'C': "00FF00", 'O': "22FF0000", 'N': "220000FF"}

traj_paths = ("/mnt/2_hdd/ter3/", "/mnt/2_hdd/gmo/", "/mnt/2_hdd/sa/")


lig_n = input("Which ligand (0: TER, 1: GMO, 2: STA): ")
traj_path = traj_paths[int(lig_n)]
starting_frame = 700 # To start calculations after relaxation

# Loading topology and trajectory files into the Universe
u = MDAnalysis.Universe(traj_path+'com-wat.top',
traj_path+'combined.nc',
traj_path+'prod1.nc',
traj_path+'prod2.nc'
)

def dist(name):
    # Measures the distance between the position of an atom in the first frame and its position in subsequent frames 
    
    l =[] # Differences in the distances of each atom is appended into this list to calculate the mean and std deviation
    
    # Choose the atom and get its position
    n = u.select_atoms(f'resname UNK and name {name}')[0]
    origin_pos = n.position


    with open(name+".dat", "w") as f:
        # Write the distances into data files to be plotted
        f.write("#"+name+"Frame\tDist\n")
        for ts in u.trajectory[starting_frame::10]:
            pos = n.position
            tot_dist = distance_array(origin_pos, pos, box=ts.dimensions)[0][0]

            l.append(tot_dist)
            f.write(f"{str(ts.frame)}\t{str(tot_dist)}\n")

    average_dist = sum(l)/len(l)
    std_dev = statistics.stdev(l) # Sample standard deviation as not all frames are included
    
    return (average_dist, std_dev)

def prep_plot(path="."):
    # Prepares the gnuplot script to plot the data files produced by dist(name) by specifying colors and writing the gnuplot script
    full_string = ""
    files = [x for x in os.listdir(path) if ".dat" in x and x[0] in ("C", "O", "N")]
    title = traj_path.split("/")[3][:3].upper()
    
    for file in files:
        name = file.strip(".dat")
        hex_color = "".join([str(ord(x)) for x in file[:3]])
        
        color = colors[name[0]]
        if name[0] == "C":
            color = "55"+rand_carbon_plot_color()

        gnuplot_string = f'"{path}/{file}" title "{name}"  lt 6 lc rgb "#{color}" lw 2 with l, \\\n'
        full_string += gnuplot_string
    
    print(full_string)
    with open("plot_all.gnu", "w") as f:
        f.write("set title '{title}' font '{{/Times:Bold boldface-newfont}}'\n".format(title=title))
        f.write("set nokey\n")
        f.write("plot ")
        f.write(full_string)
        f.write('\n pause -1 "Hit any key to continue\n"')

def write_stdev(lig_atoms):
    # Writes standard dev of atoms of a ligand to a data file
    lig_name = traj_path.split("/")[3][:3].upper()
    
    with open(lig_name+"_stdev.dat", 'w') as f:
        f.write(f"#{lig_name}_atoms\tstdev\n")
        
        for atom in lig_atoms.items():
            f.write(f"{atom[0]}\t{atom[1][1]}\n")

def rand_carbon_plot_color():
    # Generate a random shade of green for the carbon atoms in the plot
    hex_characters = "0123456789abcdef"
    sele_c = "".join(random.sample(hex_characters, 4))

    hex_color = '{}ff{}'.format(sele_c[:2], sele_c[2:])

    return hex_color

atoms = u.select_atoms('resname UNK and not (name H*)')
all_atoms = dict()
for atom in atoms:
    (avg_dist, std_dev) = dist(atom.name)
    all_atoms[atom.name] = (avg_dist, std_dev)

sorted_all_atoms = sorted(all_atoms.items(), key=lambda x:x[1][1]) # sort atoms by std deviation

list(map(lambda x:print(f"{x[0]}: {x[1][0]}, (STDEV: {x[1][1]})\n------"), sorted_all_atoms))

prep_plot() # Plots all .dat files in the current specified directory, ./ by default
write_stdev(all_atoms)