from pyteomics import mgf
filepath = "./data/9species/test_species/bacillus.10k.mgf"
spec_reader = mgf.read(filepath)
f = open("denovo_withoutPMC.txt")
denovo = list(f.readlines())

f_out = open("denovo_withoutPMC_new.txt", "w")

for idx, each in enumerate(spec_reader):
    
    c = denovo[idx].split()
    replaced_aa = c[4].replace("$", "").replace("N+0.984", "D").replace("Q+0.984", "E").replace("L","I")
    label_aa = c[1].replace("N+0.984", "D").replace("Q+0.984", "E").replace("L","I")
    c[4] = c[4].replace("$", "")
    if replaced_aa == label_aa:
        c[-1] = "correct"
    else:
        #print(replaced_aa)
        #print(c[1])
        c[-1] = "incorrect"
    final = " ".join(c)
    f_out.write(final+ "\n")
        
    

    #f_out.write(title + " " + c[1] + " " + c[4] + " " + w[4] + " " + p[4] + " " + better + " " + change + "\n")
    
    
    