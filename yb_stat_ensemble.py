import numpy as np
from pyteomics import mgf
from collections import OrderedDict

match_mass_tol = 0.1  # in Da
prefix_mass_tol = 0.5  # in Da

tools_list = [
    "Novor",
    "pNovo",
    "DeepNovo",
    "SMSNet",
    "PointNovo",
    "Casanovo"
]

_PAD = "_PAD"
_GO = "_GO"
_EOS = "_EOS"
_START_VOCAB = [_PAD, _GO, _EOS]

PAD_ID = 0
GO_ID = 1
EOS_ID = 2
vocab_reverse = [
    "A",
    "R",
    "N",
    "n",
    "D",
    "C",
    "E",
    "Q",
    "q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "m",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

vocab_reverse_nomods = [
    "a",
    "r",
    "d",
    "c",
    "e",
    "g",
    "h",
    "i",
    "l",
    "k",
    "f",
    "p",
    "s",
    "t",
    "w",
    "y",
    "v",
]

vocab_reverse = _START_VOCAB + vocab_reverse
vocab = dict([(x, y) for (y, x) in enumerate(vocab_reverse)])
vocab_size = len(vocab_reverse)

# mass value
mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949
mass_Phosphorylation = 79.96633

mass_AA = {
    "_PAD": 0.0,
    "_GO": mass_N_terminus - mass_H,
    "_EOS": mass_C_terminus + mass_H,
    "A": 71.03711,  # 0
    "R": 156.10111,  # 1
    "N": 114.04293,  # 2
    "n": 115.02695,
    "D": 115.02694,  # 3
    "C": 160.03065,  # 103.00919,  # 4
    "E": 129.04259,  # 5
    "Q": 128.05858,  # 6
    "q": 129.0426,
    "G": 57.02146,  # 7
    "H": 137.05891,  # 8
    "I": 113.08406,  # 9
    "L": 113.08406,  # 10
    "K": 128.09496,  # 11
    "M": 131.04049,  # 12
    "m": 147.0354,
    "F": 147.06841,  # 13
    "P": 97.05276,  # 14
    "S": 87.03203,  # 15
    "T": 101.04768,  # 16
    "W": 186.07931,  # 17
    "Y": 163.06333,  # 18
    "V": 99.06841,  # 19
    "x": 42.010565,
    'y': 43.005814 ,
    "z":-17.026549
}

mass_ID = [mass_AA[vocab_reverse[x]] for x in range(vocab_size)]
mass_ID_np = np.array(mass_ID, dtype=np.float32)
mass_AA_min = mass_AA["G"]  # 57.02146


def arePermutation(str1, str2):
    n1 = len(str1)
    n2 = len(str2)
    if str1 == str2:
        return False
    if n1 != n2:
        return False
    a = sorted(str1)
    str1 = " ".join(a)
    b = sorted(str2)
    str2 = " ".join(b)
    for i in range(0, n1, 1):
        if str1[i] != str2[i]:
            return False
    return True


# Function from DeepNovo to calculate correct match between two sequences
def _match_AA_novor(target, predicted):
    num_match = 0
    target_len = len(target)
    predicted_len = len(predicted)
    target_mass = [mass_ID[x] for x in target]
    target_mass_cum = np.cumsum(target_mass)
    predicted_mass = [mass_ID[x] for x in predicted]
    predicted_mass_cum = np.cumsum(predicted_mass)

    i = 0
    j = 0
    while i < target_len and j < predicted_len:
        if abs(target_mass_cum[i] - predicted_mass_cum[j]) < 0.5:
            if abs(target_mass[i] - predicted_mass[j]) < 0.1:
                num_match += 1
            i += 1
            j += 1
        elif target_mass_cum[i] < predicted_mass_cum[j]:
            i += 1
        else:
            j += 1

    return num_match


DIMENSION = 90000
BIN_SIZE = 0.1


def parse_spectra(sps):
    db = []
    for sp in sps:
        param = sp["params"]

        c = int(str(param["charge"][0])[0])

        if "seq" in param:
            pep = param["seq"]
        else:
            pep = param["title"]

        if "pepmass" in param:
            mass = param["pepmass"][0]
        else:
            mass = float(param["parent"])

        if "hcd" in param:
            try:
                hcd = param["hcd"]
                if hcd[-1] == "%":
                    hcd = float(hcd)
                elif hcd[-2:] == "eV":
                    hcd = float(hcd[:-2])
                    hcd = hcd * 500 * cr[c] / mass
                else:
                    raise Exception("Invalid type!")
            except:
                hcd = 0
        else:
            hcd = 0
        mz = sp["m/z array"]
        it = sp["intensity array"]

        db.append(
            {"pep": pep, "charge": c, "mass": mass, "mz": mz, "it": it, "nce": hcd}
        )

    return db


def readmgf(fn):
    file = open(fn, "r")
    data = mgf.read(
        file, convert_arrays=1, read_charges=False, dtype="float32", use_index=False
    )
    codes = parse_spectra(data)
    return codes

from pyteomics import mgf, mass
import pandas as pd
import os 

def fragments(peptide, types=('b', 'a', 'y')):
    """
    The function generates possible m/z for fragments of 8 ion types
    b, b2+, b-NH3, b-H2O, y, y2+, y-NH3, y-H2O
        :param
            peptide: peptide sequence as string
            types: types of fragment ions f.e. (b,y) or (a,x)
        :return
            b_ions: list of list of b_ions
            y_ions: list of list of y_ions
    """
    a_ions = []
    b_ions = []
    y_ions = []

    for i in range(1, len(peptide)):
        for ion_type in types:
            if ion_type[0] in 'b':
                b_ion_types = []
                # b
                b_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA))
                # b(2+)
                b_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=2, aa_mass=mass_AA))
                # b-H2O
                b_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='H2O', charge=1))
                # b-NH3
                b_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='NH3', charge=1))
                b_ions.append(b_ion_types)
            if ion_type[0] in 'a':
                a_ion_types = []
                # a
                a_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA))
                # a(2+)
                a_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=2, aa_mass=mass_AA))
                # a-H2O
                a_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='H2O', charge=1))
                # a-NH3
                a_ion_types.append(mass.fast_mass(
                        peptide[:i], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='NH3', charge=1))
                a_ions.append(a_ion_types)
            else:
                y_ion_types = []
                # y
                y_ion_types.append(mass.fast_mass(
                        peptide[i:], ion_type=ion_type, charge=1, aa_mass=mass_AA))
                # y(2+)
                y_ion_types.append(mass.fast_mass(
                        peptide[i:], ion_type=ion_type, charge=2, aa_mass=mass_AA))
                # y-H2O
                y_ion_types.append(mass.fast_mass(
                        peptide[i:], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='H2O', charge=1)) # if its double charged is it only half the mass
                # y-NH3
                y_ion_types.append(mass.fast_mass(
                        peptide[i:], ion_type=ion_type, charge=1, aa_mass=mass_AA) - mass.calculate_mass(formula='NH3', charge=1))
                y_ions.append(y_ion_types)
    return a_ions, b_ions, y_ions

def lcs(s1, s2):
    """Length of Longest Common Substring between two strings"""
    matrix = [["" for x in range(len(s2))] for x in range(len(s1))]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                if i == 0 or j == 0:
                    matrix[i][j] = s1[i]
                else:
                    matrix[i][j] = matrix[i - 1][j - 1] + s1[i]
            else:
                matrix[i][j] = max(matrix[i - 1][j], matrix[i][j - 1], key=len)
    cs = matrix[-1][-1]
    return len(cs)

def noise_and_fragmentIons(df, mgf_in_path):
    """Add and calculates fragment ions from a summary df
        :param
            df: summary_df
            mgf: mgf input file
        :return
            df: summary_df with additional columns 'Number of missing cleavages' and 'Median noise intensity'
            median_noise: The median noise value for all spectra in this dataset
    """  
    with mgf.read(mgf_in_path) as spectra:

        missing_cleavages = []
        cleavages_positions = []
        missing_cleavages_including_aions = []
        cleavages_positions_including_aions = []
        noise_intensities = []
        peptides=[]
        titles=[]
        denovo = open("./human/denovo_prime_new.txt")#change here
        denovo = list(denovo.readlines())
        for idx, spectrum in enumerate(spectra):
            # print(spectrum)
            #print(spectrum['params']['seq'])
            peptide = denovo[idx].split()[4]
            #print(peptide)
            peptides.append(peptide)
            title = spectrum['params']['title']
            titles.append(title)
            # peptide = peptide[::-1]
            # print(peptide)
            peptide = peptide.replace("Q+0.984",'q').replace("N+0.984",'n').replace("C+57.021","C").replace("M+15.995","m").replace("+42.011",'x').replace("+43.006","y").replace("-17.027","z")
            if type(peptide) != float:
                a_ions, b_ions, y_ions = fragments(peptide)
                # print("a_ions_length",len(a_ions))
                # print(b_ions)
                missing_cleavage = 0
                missing_cleavage_including_aions = 0
                cleavage_position = []
                cleavage_position_inlcuding_aions = []
                for pos, (a_ion, b_ion, y_ion) in enumerate(zip(a_ions, b_ions, reversed(y_ions))):
                    # with open("bacillus_ions.txt","a") as f:
                    #     f.write(str(spectrum["params"]["title"])+"\n")
                    # print("one_a:",a_ion)
                    a_ions_present = [True for i in a_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    b_ions_present = [True for i in b_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    y_ions_present = [True for i in y_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    if b_ions_present or y_ions_present:
                        cleavage_position.append(pos+1)
                    else:
                        missing_cleavage +=1
                    
                    if b_ions_present or y_ions_present or a_ions_present:
                        cleavage_position_inlcuding_aions.append(pos+1)
                    else:
                        missing_cleavage_including_aions += 1
                missing_cleavages.append(int(missing_cleavage))
                cleavages_positions.append(cleavage_position)

                missing_cleavages_including_aions.append(int(missing_cleavage_including_aions))
                cleavages_positions_including_aions.append(cleavage_position_inlcuding_aions)


                spectrum_intensity = spectrum['intensity array'] 
                total_b = 0
                total_y = 0
                for pos, (b_ion, y_ion) in enumerate(zip(b_ions, reversed(y_ions))):
                    b = [np.nonzero(np.isclose(spectrum['m/z array'], i, atol=match_mass_tol)) for i in b_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    y = [np.nonzero(np.isclose(spectrum['m/z array'], i, atol=match_mass_tol)) for i in y_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    
                    for elem in b:
                        
                        for index in elem:
                            total_b += len(index.tolist())
                            #with open("case_bions.txt","a") as f:
                                #f.write("B_ions:"+ str(spectrum['m/z array'][index]) + "\n")
                        spectrum_intensity[elem] = np.nan
                    for elem in y:
                        
                        for index in elem:
                            total_y += len(index.tolist())
                            #with open("case_yions.txt","a") as f:
                                #f.write("Y_ions:"+ str(spectrum['m/z array'][index]) + "\n")
                        spectrum_intensity[elem] = np.nan
                #print(total_y)
                #print(total_b)
                with open("./human/prime_by_scores2.txt", "a") as f:
                    f.write(str(total_y+total_b) + "\n")
                noise_intensities = noise_intensities + list(spectrum_intensity)
            else:
                missing_cleavages.append(np.nan)
                cleavages_positions.append(np.nan)
                missing_cleavages_including_aions.append(np.nan)
                cleavages_positions_including_aions.append(np.nan)
        df["PepLabel"] = peptides
        df["titles"] = titles
        df["Number of missing cleavages"] = missing_cleavages
        df["Position of present cleavages"] = cleavages_positions
        df["Number of missing cleavages (including a-ions)"] = missing_cleavages_including_aions
        df["Position of present cleavages (including a-ions)"] = cleavages_positions_including_aions
        median_noise = np.nanmedian(noise_intensities)
        median_missing_cleavages = np.nanmedian(missing_cleavages)
        # logger.info(f"Median noise intensity: {median_noise}")
        # logger.info(f"Median number of missing cleavages: {median_missing_cleavages}")
        return df, median_noise
    
def noise_factor(df, mgf_in_path, median_noise):
    """Add and calculates noise factor from a summary df

        :param
            df: summary_df
            mgf: mgf input file
        :return
            df: dataframe with additional columns 'Noise factor', 'Number of noise peaks', 'Number of fragment peaks'
    """
    noise_factor_list = [] 
    noise_total_amount = []
    fragments_total_amount = []

    with mgf.read(mgf_in_path) as spectra:
        for spectrum in spectra:
            peptide = spectrum['params']['seq']
            peptide = peptide.replace("Q+0.984",'q').replace("N+0.984",'n').replace("C+57.021","C").replace("M+15.995","m").replace("+42.011",'')
            if type(peptide) != float:
                amount_noise_peaks = 0
                amount_fragment_peaks = 0
                a_ions, b_ions, y_ions = fragments(peptide)
                spectrum_intensity = spectrum['intensity array'] 
                for pos, (b_ion, y_ion) in enumerate(zip(b_ions, reversed(y_ions))):
                    b = [np.nonzero(np.isclose(spectrum['m/z array'], i, atol=match_mass_tol)) for i in b_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    y = [np.nonzero(np.isclose(spectrum['m/z array'], i, atol=match_mass_tol)) for i in y_ion if np.isclose(spectrum['m/z array'], i, atol=match_mass_tol).any()]
                    for elem in b:
                        spectrum_intensity[elem] = np.nan
                        amount_fragment_peaks += 1
                    for elem in y:
                        spectrum_intensity[elem] = np.nan
                        amount_fragment_peaks += 1
                
                for peak in spectrum_intensity:
                    if peak >= median_noise:
                        amount_noise_peaks += 1

                noise_total_amount.append(amount_noise_peaks)
                fragments_total_amount.append(amount_fragment_peaks)
                if amount_fragment_peaks > 0:
                    noise_val = amount_noise_peaks/amount_fragment_peaks
                else:
                    noise_val = amount_noise_peaks/1
                noise_factor_list.append(noise_val)
            else:
                noise_factor_list.append(np.nan)
                fragments_total_amount.append(np.nan)
                noise_total_amount.append(np.nan)
        df["Noise factor"] = noise_factor_list
        df["Number of noise peaks"] = noise_total_amount
        df["Number of fragment peaks"] = fragments_total_amount

        total_peaks = np.nansum(noise_total_amount) + np.nansum(fragments_total_amount)
        noise_percent = np.nansum(noise_total_amount) / total_peaks
        # logger.info(f"{noise_percent*100}% of all peaks are noise")
        # logger.info(f"Average noise factor: {np.nanmean(noise_factor_list)}")
        return df

# species = ["bacillus","mouse","human","ricebean","yeast","clambacteria","mmazei","tomato","honeybee"]

species = ["bacillus"]

for specie in species:
    df = pd.DataFrame(data={})
    dataPath = "../data/9species/test_species/human.10k.mgf"
    df, median_noise = noise_and_fragmentIons(df,dataPath)
    # missing = df["Number of missing cleavages"]
    # df = noise_factor(df,dataPath,median_noise)
    # df.to_csv(specie + "_results.csv",index=True)
