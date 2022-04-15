import os

def read_text_file(file_path):
    with open(file_path,"r") as f:
        # print(f.read())
        lat_const, total_energy = "",""

        flag_exclamation = False
        
        for l in f.readlines():
            if "lattice parameter (alat)  =" in l:
                lat_const = float(l.lstrip("lattice parameter (alat)  = ").split()[0])
            
            if "!" in l:
                total_energy = l.lstrip("!    total energy              = ").rstrip(" Ry\n") + "\n"
                flag_exclamation =True
    if flag_exclamation:
        return lat_const, total_energy
    else: 
        return False


if __name__ == '__main__':
    # edit path to the dir with results from scf 
    path = "/Users/adamsmith/Documents/GitHub/PH30036_AdamSmith/FYP/Ta2CO2/ConvergenceTest/lat_const/lat_const_CT"
    outputfile_end = ".out"
    # change dir to the path above
    os.chdir(path)

    # create a new file to put data in!
    new_f = open("data.txt", "w")
    new_f.write("Results for the lattice constant convegence in Ta2CO2\n\nlattice constant (bohr)           total energy  (Ry)\n")
    

    data = []
    # iterate through all files in file
    for file in os.listdir():
        # check format of the file
        if file.endswith(outputfile_end):
            # n the value entered in the scf input for generating the unit cell calculated over
            n = file.split(".")[2]
            file_path = path + "/" + file
            output = read_text_file(file_path)
            # check to make sure the .out had a total energy calculated
            if output !=False:
                lat_const, totalE = output
                data.append([lat_const,totalE])
                # new_f.write(str(n + "       " + lat_const + "           " + totalE))
    data.sort(key=lambda x:float(x[0]))  
    for d in data:
        new_f.write(str(str(d[0]) + "       " + str(d[1])))  
    new_f.close()