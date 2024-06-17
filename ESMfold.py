import torch
import esm
import os
from Bio import SeqIO

import mutations

print("This program will automatically process the mutation in interest and return predicted PDB file using ESMfold")
print("To start with, process the fasta file for your mutation tagets")
fasta_continue = "y"
while fasta_continue == "y":
    mutations.main()
    fasta_continue = input('DO you need to process another fasta file? [y/n]')

print("Now we will help you load the ESMfold module in your bridge directory\n")
print("------------------------------Notice--------------------------------------")
print("Please make sure to run the model with an interact node and have your pwd located in Ocean, otherwise program will be killed due to insufficient space")
print("Your current working directory is:", os.getcwd(), "\n")
pwd_continue = input("Do you want to change your working directory? [y/n]")
if pwd_continue == "y":
    newpwd = ""
    try:
        newpwd = input("Please enter your new working directory:\n")
        os.chdir(newpwd)
    except FileNotFoundError:
        print("The directory does not exist")
    except PermissionError:
        print("You do not have permissions to change to this directory")
    except Exception as e:
        print(f"An error occurred: {e}")
print("Loading ESMfold Module, Please wait...\n")
os.environ["TORCH_HOME"] = os.getcwd() + "/torch_cache"
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()

print("Now, we will help you to run the PDB model prediction.\n")

def fasta_read():
    while True:
        fasta_path = input("Please input the path for the fasta file that you want to predict:\n")
        try:
            records = SeqIO.parse(fasta_path, "fasta")
            for record in records:
                sequence = str(record.seq)
            break
        except Exception:
            print("ERROR: Destinated file not found!\n Please check if the path is correct, then re-enter the path again.\n")
            fasta_read()
    return sequence

esm_continue = "y"
while esm_continue == "y":
    sequence = fasta_read()
    with torch.no_grad():
        print(sequence)
        output = model.infer_pdb(sequence)
    esm_save = input("Please enter the file name you want to save as, filename should end with .pdb:\n")
    with open("ESM_predictions/"+ esm_save, "w") as f:
        f.write(output)
    print("Your file has been saved as:", esm_save, "in the subdirectory /ESM_predictions\n\n")
    esm_continue = input('Do you need to predict another fasta file? [y/n]')
print("-------------------------------------------------------------------------------------------------------------------------------------------------------")
print("Due to module compatibility issue, the analysis of pdb files can not be completed on bridges. ")
print("We offer instruction for RMSD and secondary structure assignment comparison for your PDB result. Please refer to the user manual for how to set up the environment locally and run the PDB_comparison.py ")
print("Thanks for using!")
print("-------------------------------------------------------------------------------------------------------------------------------------------------------")



