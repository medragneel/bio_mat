from math import log2
import argparse

parser=argparse.ArgumentParser(description='Bio_mat')
parser.add_argument('-f','--file',help='read from file')
parser.add_argument('-s','--seqs',nargs="*",help='read seqs from inputs')
args= parser.parse_args()
print(args.seqs)
# creer une matrice de nrows (lines) et ncols (columns)
def create_mat(nrows,ncols):
    mat = [[0.0] * ncols for _ in range(0,nrows)]
    return mat
def print_mat(mat,alphabet="ACGT"):
    for i in range(len(mat)):
        print(alphabet[i],end="\t")
        for j in range(len(mat[i])):
            print(str(mat[i][j]).rstrip("0").rstrip("."),end="\t")
        print()


    
def mat_counts(seqs,alphabet="ACGT"):
    counts = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            counts[idx][i] += 1
    return counts


def fq_mat(seqs,alphabet="ACGT"):
    c = mat_counts(seqs)
    f = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            f[idx][i] = float("{0:.2f}".format(c[idx][i] / len(seqs)))
    return f

def pwm(seqs,alphabet="ACGT"):
    f  = fq_mat(seqs)
    p = create_mat(len(alphabet),len(seqs[0]))
    for s in seqs:
        for i in range(len(seqs[0])):
            idx = alphabet.index(s[i])
            p[idx][i] =  float("{0:.2f}".format(log2(f[idx][i] / 0.25)))
    return p 
    

def get_motif(seqs,c,alphabet="ACGT"):
    m = ""
    for j in range(len(seqs[0])):
        max_c = c[0][j]
        max_ci = 0
        for i in range(1,len(alphabet)):
            if c[i][j] > max_c:
                max_c = c[i][j]
                max_ci =  i
        m += alphabet[max_ci] 
    return m

def read_seq_form_file(file):
    with open(file,"r") as f:
        lines = [line.strip("\n") for line in f.readlines() ]
    return lines



seqs2=[
        "GCCGGAAGTG",
        "ACCGGAAGCA",
        "GCCGGATGTA",
        "ACCGGAAGCT",
        "ACCGGATATA",
        "CCCGGAAGTG",
        "ACAGGAAGTC",
        "GCCGGATGCA",
        "TCCGGAAGTA",
        "ACAGGAAGCG",
        "ACAGGATATG",
        "TCCGGAAACC",
        "ACAGGATATC",
        "CAAGGACGAC",
        "TCTGGACCCT",
        ]

seqs = [
        "AAAGTT",
        "CACGTG",
        "TTGGGT",
        "GACCGT",
        "AACCAT",
        "AAACCT",
        "GAACCT"
        ]
def test():
    if args.file:
        s = read_seq_form_file(args.file)
        c = mat_counts(s)
        f = fq_mat(s)
        p = pwm(s)
        print("matrice de contage")
        print_mat(c)
        print("matrice de frequence")
        print_mat(f)
        print("pwm")
        print_mat(p)
        print("motif")
        m= get_motif(s,p)
        print(m)
    elif args.seqs:
        c = mat_counts(args.seqs)
        f = fq_mat(args.seqs)
        p = pwm(args.seqs)
        print("matrice de contage")
        print_mat(c)
        print("matrice de frequence")
        print_mat(f)
        print("pwm")
        print_mat(p)
        print("motif")
        m= get_motif(args.seqs,p)
        print(m)

if __name__ == "__main__":
    test()
