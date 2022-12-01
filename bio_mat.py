from math import log2

def log(x,base):
    result = ln(x)/ln(base)
    return result

def ln(x):
    n = 1000.0
    return n * ((x ** (1/n)) - 1)

def create_mat(nrows,ncols):
    mat = [[0.0] * ncols for _ in range(0,nrows)]
    return mat
def print_mat(mat,alphabet="ACGT"):
    for i in range(len(mat)):
        print(alphabet[i],end="\t")
        for j in range(len(mat[i])):
            print(mat[i][j],end="\t")
        print()

class Bio_mat:
    def __init__(self,seqs,alphabet="ACGT"):
        self.alphabet = alphabet
        self.seqs = seqs
    def mat_counts(self):
        counts = create_mat(len(self.alphabet),len(self.seqs[0]))
        for s in self.seqs:
            for i in range(len(self.seqs[0])):
                idx = self.alphabet.index(s[i])
                counts[idx][i] += 1
        return counts

    def fq_mat(self):
        c = self.mat_counts()
        f = create_mat(len(self.alphabet),len(self.seqs[0]))
        for s in self.seqs:
            for i in range(len(self.seqs[0])):
                idx = self.alphabet.index(s[i])
                f[idx][i] = float("{0:.2f}".format(c[idx][i] / len(self.seqs)))
        return f
    def pwm(self):
        f  = self.fq_mat()
        p = create_mat(len(self.alphabet),len(self.seqs[0]))
        for s in self.seqs:
            for i in range(len(self.seqs[0])):
                idx = self.alphabet.index(s[i])
                p[idx][i] =  float("{0:.2f}".format(log2(f[idx][i] / 0.25)))
        return p 

    def get_motif(self,c):
        m = ""
        for j in range(len(self.seqs[0])):
            max_c = c[0][j]
            max_ci = 0
            for i in range(len(self.alphabet)):
                if c[i][j] > max_c:
                    max_c = c[i][j]
                    max_ci =  i
            m += self.alphabet[max_ci] 
        return m
        


def test():
    
    seqs=[
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
    bm = Bio_mat(seqs)
    c = bm.mat_counts()
    f = bm.fq_mat()
    p = bm.pwm()
    print("matrice de contage")
    print_mat(c)
    print("matrice de frequence")
    print_mat(f)
    print("matrice PWM")
    print_mat(p)
    print("motif")
    print(bm.get_motif(p))

    
if __name__ == "__main__":
    test()
