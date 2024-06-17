import Bio.PDB
from Bio.PDB import DSSP
import csv
#conda install -c salilab dssp

def RMSD_process(sample_model, ref_model, sample_file, ref_file):

  # Make a list of the atoms (in the structures) you wish to align.
  ref_atoms = []
  sample_atoms = []

  for sample_chain in sample_model:
    for sample_res in sample_chain:
      if(sample_res.has_id("CA")):
          sample_atoms.append(sample_res["CA"])
      else:
        print("No carbon found in", sample_res, "Please check your PDB file")
        return None
  for ref_chain in ref_model:
    for ref_res in ref_chain:
        if(ref_res.has_id("CA")):
            ref_atoms.append(ref_res["CA"])
        else:
          print("No carbon found in", sample_res, "Please check your PDB file")
          return None
    if (len(ref_atoms) != len(sample_atoms)):
        print("The length of your two sequence are not the same. Please check your model")
        return None
    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Print RMSD:
    print ("RMSD for alignment: Reference Sample:", ref_file, "Align Sample:", sample_file,"\n", super_imposer.rms)

def Secondary_Structure_Process(sample_model, ref_model, sample_file, ref_file):
    structure_assignment = []

    # Run DSSP analysis - replace 'dssp' with the path to the DSSP executable if not in PATH
    ref_dssp = DSSP(ref_model, ref_file, dssp='mkdssp')
    sample_dssp = DSSP(sample_model, sample_file, dssp='mkdssp')
    
    # Print DSSP results for each residue
    print("Analyzing Secondary Structure Assignment for PDB file", ref_file)

    #Store assignment of each model in list for output
    for residue in ref_dssp:
        index, amino_acid, secondary_structure = residue[0], residue[1], residue[2]
        assignment = {}
        assignment["index"] = index
        assignment["Reference_amino_acid"] = amino_acid
        assignment["Reference_secondary_structure"] = secondary_structure
        structure_assignment.append(assignment)
    print("Analyzing Secondary Structure Assignment for PDB file", sample_file)
    for residue in sample_dssp:
        index, amino_acid, secondary_structure = residue[0], residue[1], residue[2]
        assignment = structure_assignment[index-1]
        assignment["Sample_amino_acid"] = amino_acid
        assignment["Sample_secondary_structure"] = secondary_structure
    print("Finish Assigning structures for both files.")
    
    #Output delivery
    option = input("How do you wish to check the assignment?\n[1]Save to local csv file\n[2]Print out in terminal\n[3]Do nothing")
    if option == "1":
        filename = input("input the file name you want to save as, file name should end with .csv:")
        filename = filename
        # Open the file in write mode
        with open(filename, 'w', newline='') as csvfile:
            # Identify the fieldnames from the keys of the first dictionary
            fieldnames = structure_assignment[0].keys()
            # Create a DictWriter object with the fieldnames
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            # Write the header (the field names)
            writer.writeheader()
            for data in structure_assignment:
                writer.writerow(data)
        csvfile.close()
    elif option == "2":
        for row in structure_assignment:
            print(row)
# Example usage


#Global Variables
def PDB_process():
    print("This program will help you to analyze the RMSD and the secondary structure assignment of your referene and mutated sample structure")
    print("Note: This program CAN NOT process pdb models with different length of sequence! Please make sure your input model are for same length, or the program will terminate when it finds the mismatch error")
    ref_file = input("Please input the path for your reference model:")
    sample_file = input("Please input the path for your sample model (model to be aligned):")

    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)

    # Get the structures
    ref_structure = pdb_parser.get_structure("reference", ref_file)
    sample_structure = pdb_parser.get_structure("samle", sample_file)

    # Use the first model in the pdb-files for alignment
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]

    RMSD_process(sample_model, ref_model, sample_file, ref_file)
    Secondary_Structure_Process(sample_model, ref_model, sample_file, ref_file)

PDB_process()    