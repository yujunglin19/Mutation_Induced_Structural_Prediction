import os
def user_input():
    input_fasta = input("Enter fasta file path:\n")
    try:
        with open(input_fasta, 'r') as f:
            seq_id = ""
            aa_seq = ""
            for line in f:
                if '>' in line:
                    tmp = list(line.strip().split('|'))
                    seq_id = tmp[0][1:]
                else:
                    aa_seq += line.strip()
        path = os.path.dirname(input_fasta)
        #print(path)
    except Exception:
        print("ERROR: Destinated file not found!\n Please check if the path is correct, then re-enter the path again.")
        user_input()

    return seq_id, aa_seq, path

def input_mutation():
    input_mut = input("Enter mutations: (if inputing multiple mutations, separate each by space)\n")
    mutations = list(input_mut.strip().split(' '))
    return mutations

def mutation(mutations, aa_seq):
    performed_mut = []
    for m in mutations:
        if m[0].isdigit():
            pos = int(m[:-1])
            s = m[-1]
        else:
            o = m[0]
            pos = int(m[1:-1])
            s = m[-1]
        if pos > len(aa_seq):
            print(f"Mutation %s skipped. Input mutation not in the range of input protein sequence." % m)
            print(' ')
        else:
            s = s.upper()
            og_aa = aa_seq[pos-1]
            aa_seq = aa_seq[:pos-1] + s + aa_seq[pos:]
            if m[0].isdigit():
                m = og_aa + m
            else:
                m = og_aa.upper() + m[1:]
            performed_mut.append(m)
    return aa_seq, performed_mut

def output_file (path, aa_seq, seq_id, mutations):
    user_input = input("Save sequence to fasta file: [y/n]\n")
    if user_input.strip() == 'n':
        return
    elif user_input.strip() == 'y':
        file_name = seq_id + '_' + '_'.join(mutations) + '.fasta'
        print(file_name)
        save_path = os.path.join(path, file_name)
        print(save_path)
        first_line = '>mut_' + seq_id + '_' + '_'.join(mutations)
        print(first_line)
        with open(save_path, "w") as f:
            f.write(first_line + '\n')
            f.write(aa_seq)
        print('Output fasta file saved to ' + save_path)
    else:
        print("User input not recognized. Please re-enter the response.\n")
        output_file(aa_seq, seq_id, mutations)

def mutation_option(seq_id, aa_seq, path):
    ans = input("Do you want to perform mutations? [y/n]\n")
    if ans == 'y':
        mutations = input_mutation()
        mut_aa_seq, performed_mut = mutation(mutations, aa_seq)
        print(f"Protein ID: %s" % seq_id)
        print(' ')
        print(f"Original protein sequence: \n%s" % aa_seq)
        print(' ')
        print("Input mutations:\n")
        print(' '.join(mutations))
        print(' ')
        print("Mutation performed:\n")
        print(' '.join(performed_mut))
        print(' ')
        print("Mutated protein sequence:\n")
        print(mut_aa_seq)
        output_file(path, mut_aa_seq, seq_id, performed_mut)
    elif ans == 'n':
        print('Mutation skipped. Continue to structural prediction.')
    else:
        print("User input not recognized. Please re-enter the response.\n")
        mutation_option(seq_id, aa_seq, path)

def main():
    seq_id, aa_seq, path = user_input()
    mutation_option(seq_id, aa_seq, path)



main()
