from pyteomics import mgf
import re
species = "human"
def calc(sequence):
    sequence = sequence.replace("I", "L")
    sequence = re.split(r"(?<=.)(?=[A-Z])", sequence)
    total = 0
    for each in sequence:
        
        # total += aa2mas[each]
        
        try:
            total += residues[each]
        except:
            h1 = each.count("+42.011")
            h2 = each.count("+43.006")
            h3 = each.count("-17.027")
            total += h1 * 42.010565 + h2 * 43.005814  + h3 * -17.026549
            each = each.replace("+42.011", "")
            each = each.replace("+43.006", "")
            each = each.replace("-17.027", "")
            if each:
                total += residues[each]
                    
    return total 
#model_new is the washed model output with correct/incorrect label
prime = open("{}/denovo_prime_new.txt".format(species))
prime = list(prime.readlines())
casa = open("./denovo_casa_new.txt")
casa = list(casa.readlines())
clip = open("{}/denovo_clip_new.txt".format(species))
clip = list(clip.readlines())
without =  open("./denovo_withoutPMC_new.txt")
without = list(without.readlines())
reverse = open("{}/denovo_reverse_clip_new.txt".format(species))
reverse = list(reverse.readlines())

#scores_xx is the confidence score by model
scores = open("{}/score_reverse.txt".format(species))
scores_r = list(scores.readlines())
scores = open("{}/score_clip.txt".format(species))
scores_c = list(scores.readlines())
scores = open("{}/score_prime.txt".format(species))
scores_p = list(scores.readlines())

#xx_by_scores2 is the new yb scores we calculated 
yb_r = open("{}/reverse_by_scores2.txt".format(species))
yb_r = list(yb_r.readlines())
yb_c = open("{}/clip_by_scores2.txt".format(species))
yb_c = list(yb_c.readlines())
yb_p = open("{}/prime_by_scores2.txt".format(species))
yb_p = list(yb_p.readlines())
#xx_by_score_abx is new scores based on all 6 ions
# yb_r = open("reverse_by_scores_abx.txt")
# yb_r = list(yb_r.readlines())
# yb_c = open("clip_by_scores_abx.txt")
# yb_c = list(yb_c.readlines())
# yb_p = open("prime_by_scores_abx.txt")
# yb_p = list(yb_p.readlines())

#clipscores:
f = open("prime_clipscores.txt")
prime_scores = f.readlines()

f = open("reverse_clipscores.txt")
reverse_scores = f.readlines()
f = open("clip_clipscores.txt")
clip_scores = f.readlines()


import yaml

with open("config.yaml") as f:
    data = yaml.load(f,Loader=yaml.FullLoader)
if species == ".":
    filepath = "../data/9species/test_species/bacillus.10k.mgf"
else:
    filepath = '../data/9species/test_species/{}.10k.mgf'.format(species)
spec_reader = mgf.read(filepath)
residues = data["residues"]

f_out = open("{}/denovo_ensemble.txt".format(species), "w")

non_prime = 0
##prepare data for LR:
data = []
for idx, spectrum in enumerate(spec_reader):
    try:
        p_cs = float(prime_scores[idx])
        r_cs = float(reverse_scores[idx])
        c_cs = float(clip_scores[idx])
    except:
        p_cs = 0
        r_cs = 0
        c_cs = 0
        
    
    
    #ca = casa[idx].split()
    p = prime[idx].split()
    r = reverse[idx].split()
    cl = clip[idx].split()
    #w = without[idx].split()
    
    
    #print(spectrum["params"]["charge"])
    #title = spectrum['params']['title']
    # print(title)
    m_z = spectrum["params"]["pepmass"][0]
    spectrum["params"]["pepmass"]=(m_z,'')
    # spectrum["params"]["pepmass"][1] = ''
    pep = spectrum["params"]["seq"]
    charge = int(spectrum["params"]["charge"][0])
    #pep = pep.replace("C(Carbamidomethyl)", "C+57.021").replace('C+57.021','C').replace('C','C+57.021').replace("M(Oxidation)","M+15.995").replace("M(ox)","M+15.995")
    spectrum["params"]["seq"] = pep
    #charge = int(truemass / float(m_z) + 1)
    specMass =  (m_z - 1.007276) * charge - 18.01
    #mass_ca = calc (ca[4])
    mass_p = calc(p[4])
    mass_r = calc(r[4])
    mass_cl = calc(cl[4])
    #mass_w = calc(w[4])
    offset = 0.65 # offset for prime cf. 0.65
    cf_scaler = 1 #sclaer for weight of confidence score 1
    scaler = 0.01 #scaler for weight of yb ions count 0.01
    scaler_cs = 1 #scaler for weight of clip scores 1
    confidence_r = float(scores_r[idx].split()[1]) * cf_scaler 
    confidence_c = float(scores_c[idx].split()[1]) * cf_scaler
    confidence_p = float(scores_p[idx].split()[1]) * cf_scaler
    new_confidence_p = ( confidence_p + offset ) * cf_scaler
    ion_r = int(yb_r[idx])
    ion_c = int(yb_c[idx])
    ion_p = int(yb_p[idx])
    

    ion_diff_tol = 30
    
    max_score = max(confidence_c + scaler * ion_c + scaler_cs * c_cs , confidence_r + scaler * ion_r + scaler_cs * r_cs, new_confidence_p + scaler * ion_p + scaler_cs * p_cs)
    sorted_ions = sorted([ion_c, ion_r, ion_p])
    ther = 0.03
    
    if abs(mass_r - specMass) < ther and abs(mass_cl - specMass) < ther and abs(mass_p - specMass) < ther:
        if sorted_ions[-1] - sorted_ions[-2] > ion_diff_tol: 
            if sorted_ions[-1] == ion_c:
                f_out.write(clip[idx])
                print(1)
                non_prime +=1
            elif sorted_ions[-1] == ion_r:
                f_out.write(reverse[idx])
                print(2)
                non_prime +=1
            else:
                f_out.write(prime[idx])
                #print(3)
        else:
                
            if max_score == confidence_c + scaler * ion_c + scaler_cs * c_cs :
                f_out.write(clip[idx])
                non_prime+=1
                print(1)
            elif max_score ==confidence_r + scaler * ion_r + scaler_cs * r_cs:
                f_out.write(reverse[idx])
                non_prime+=1
                print(2)
            else:
                f_out.write(prime[idx])
                
                #print(3)
    elif abs(mass_r - specMass) < ther and abs(mass_cl - specMass) < ther:
        
        
            
        if confidence_r + scaler * ion_r + scaler_cs * r_cs > confidence_c + scaler * ion_c + scaler_cs * c_cs:
            f_out.write(reverse[idx])
            non_prime+=1
            print(2)
        else:
            f_out.write(clip[idx])
            non_prime+=1
            print(1)
    elif abs(mass_r - specMass) < ther and abs(mass_p - specMass) < ther:
        if ion_r > ion_p + ion_diff_tol:
            f_out.write(reverse[idx])
            non_prime+=1
            print(2)
        elif ion_p > ion_r + ion_diff_tol:
            
            f_out.write(prime[idx])
            #print(3)
        else:
            if confidence_r + scaler * ion_r + scaler_cs * r_cs > new_confidence_p + scaler * ion_p + scaler_cs * p_cs:
                f_out.write(reverse[idx])
                non_prime+=1
                print(2)
            else:
                f_out.write(prime[idx])
                #print(3)
    elif abs(mass_p - specMass) < ther and abs(mass_cl - specMass) < ther:
        if ion_c > ion_p + ion_diff_tol:
            f_out.write(clip[idx])
            print(1)
            non_prime+=1
        elif ion_p > ion_c + ion_diff_tol:
            f_out.write(prime[idx])
            #print(3)
        else:
            if new_confidence_p + scaler * ion_p + scaler_cs * p_cs > confidence_c+ scaler * ion_c + scaler_cs * c_cs:
                f_out.write(prime[idx])
                #print(3)
            else:
                f_out.write(clip[idx])
                non_prime+=1
                print(1)
    else:
        f_out.write(prime[idx])
        #print(3)
    
print(non_prime)
   
            
        
        
    
    
    
    
    #------------------------
    #ensemble algorithm
    # if abs(mass_r - specMass) < 0.02 and abs(mass_cl - specMass) < 0.02 and abs(mass_p - specMass) < 0.02:
    #     if new_confidence_p > 0.96 or confidence_r > 0.96 or confidence_c > 0.96:
    #         max_score = max(confidence_c, confidence_r, new_confidence_p)
    #         if max_score == confidence_c:
    #             f_out.write(clip[idx])
    #         elif max_score ==confidence_r:
    #             f_out.write(reverse[idx])
    #         else:
    #             f_out.write(prime[idx])
    #     else:
            
                
                
            
            
                
        
        
    # elif abs(mass_r - specMass) < 0.02 and abs(mass_cl - specMass) < 0.02:
    #     if confidence_r > 0.95 or confidence_c > 0.95:
    #         if confidence_r>confidence_c:
    #             f_out.write(reverse[idx])
    #         else:
    #             f_out.write(clip[idx])
    #     else: # if both confidence < 0.99
    #         f_out.write(prime[idx])
    # elif abs(mass_r - specMass) < 0.01:
    #     if confidence_r > 0.985:
    #         f_out.write(reverse[idx])
    #     else:
    #         f_out.write(prime[idx])
    # elif abs(mass_cl - specMass) < 0.01:
    #     if confidence_c > 0.985:
    #         f_out.write(clip[idx])
    #     else:
    #         f_out.write(prime[idx])
    # else:
    #     f_out.write(prime[idx])
        
        
            
        
            
            
        
            
            
                
            
        
    # if abs(mass_r - specMass) < 0.01 and confidence_r > 0.995:
    #     f_out.write(reverse[idx])
        
    #     print(1)
    #     # print(without[idx])
    #     # print(mass_r)
    #     # print(specMass)
    # elif abs(mass_cl - specMass) < 0.01 and confidence_c>0.995:
    #     f_out.write(clip[idx])
    #     print(2)
    # else:
    #     f_out.write(prime[idx])
    #     print(3)
        
    





    
    
    