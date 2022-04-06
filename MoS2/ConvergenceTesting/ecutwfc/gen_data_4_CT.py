import os

def read_text_file(file_path):
    with open(file_path,"r") as f:
        # print(f.read())
        num_of_kpts, total_energy = "",""

        flag_exclamation = False
        
        for l in f.readlines():
            if "number of k points= " in l:
                num_of_kpts = int(l.lstrip("number of k points= ").split()[0])
            
            if "!" in l:
                total_energy = l.lstrip("!    total energy              = ").rstrip(" Ry\n") + "\n"
                flag_exclamation =True
    if flag_exclamation:
        return num_of_kpts, total_energy
    else: 
        return False


if __name__ == '__main__':
    # edit path to the dir with results from scf 
    path = "/Users/adamsmith/Documents/GitHub/PH30036_AdamSmith/FYP/MoS2/ConvergenceTesting/ecutwfc/ecutwfc_CT"
    outputfile_end = ".out"
    # change dir to the path above
    os.chdir(path)

    # create a new file to put data in!
    new_f = open("data.txt", "w")
    
    

    data = []
    # iterate through all files in file
    for file in os.listdir():
        # check format of the file
        if file.endswith(outputfile_end):
            # n the value entered in the scf input for generating the unit cell calculated over
            ecutwfc = file.split(".")[2]
            file_path = path + "/" + file
            output = read_text_file(file_path)
            # check to make sure the .out had a total energy calculated
            if output !=False:
                kpts,totalE = output
                data.append([ecutwfc,totalE])
                # new_f.write(str(n + "       " + numofkpts + "           " + totalE))
    data.sort(key=lambda x:int(x[0]))  
    new_f.write("Results for the ecutwfc convegence in Ta2C\n"+ "For number of kpts = "+ str(kpts)  + "\necutwfc          total energy / Ry\n")
    for d in data:
        new_f.write(str(str(d[0]) + "            " + str(d[1])))  
    new_f.close()