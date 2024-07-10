import os

for file in os.listdir():
    if file[-5:] == ".xlsx":
        print(file[:-5])
        if file[-6] == "C":
            new_name_C = file[:-6] + "1" + ".xlsx"
            os.rename(file, new_name_C)
        elif file[-6] == "N":
            new_name_N = file[:-6] + "2" + ".xlsx"
            os.rename(file, new_name_N)
