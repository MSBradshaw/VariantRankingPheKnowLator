import networkx as nx
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_coms(file):
    coms = {}
    for line in open(file, 'r'):
        row = line.strip().split('\t')
        for com in row[1:]:
            if com in coms:
                coms[com].append(row[0])
            else:
                coms[com] = [row[0]]
    return coms


m_coms = get_coms('CommunityDetection/markov.txt')
l_coms = get_coms('CommunityDetection/louvain.txt')
w_coms = get_coms('CommunityDetection/walktrap.txt')

cftr = {'CFTR', 'HP:0011949', 'HP:0001733', 'HP:0002202', 'HP:0030830', 'HP:0000006', 'HP:0011962', 'HP:0002105',
        'HP:0000007', 'HP:0004313', 'HP:0002027', 'HP:0002721', 'HP:0005952', 'HP:0006528', 'HP:0006538',
        'HP:0002110', 'HP:0030877', 'HP:0002613', 'HP:0011961', 'HP:0001738', 'HP:0030247', 'HP:0002783',
        'HP:0000118', 'HP:0000837', 'HP:0005376', 'HP:0001508', 'HP:0001658', 'HP:0003251', 'HP:0000819',
        'HP:0012873', 'HP:0100749', 'HP:0002099', 'HP:0005213', 'HP:0002024', 'HP:0004401', 'HP:0000798',
        'HP:0001425', 'HP:0002205', 'HP:0100027', 'HP:0012379', 'HP:0001945', 'HP:0004469', 'HP:0011947',
        'HP:0002240', 'HP:0002150', 'HP:0011227', 'HP:0001974', 'HP:0004326', 'HP:0002570', 'HP:0001648',
        'HP:0002035'}

htt = {'HTT', "HP:0000738", "HP:0000752", "HP:0003324", "HP:0002119", "HP:0011448", "HP:0006855", "HP:0031473",
            "HP:0001336", "HP:0002529", "HP:0002376", "HP:0002340", "HP:0000716", "HP:0003487", "HP:0045082",
            "HP:0000006", "HP:0002073", "HP:0002375", "HP:0001263", "HP:0000007", "HP:0000496", "HP:0002317",
            "HP:0002591", "HP:0002066", "HP:0002312", "HP:0000718", "HP:0002355", "HP:0002300", "HP:0001262",
            "HP:0008936", "HP:0200136", "HP:0000708", "HP:0002059", "HP:0002808", "HP:0003763", "HP:0000741",
            "HP:0001824", "HP:0000734", "HP:0010864", "HP:0001773", "HP:0001288", "HP:0002540", "HP:0002650",
            "HP:0007256", "HP:0002169", "HP:0200147", "HP:0002141", "HP:0000746", "HP:0000737", "HP:0031845",
            "HP:0002354"}

aa = {'SBDS', 'NBN', 'PRF1', 'IFNG'}
tof = {'NKX2-5', 'GATA4', 'ZFPM2', 'GATA6', 'JAG1', 'TBX1'}


def get_non_zero_com_counts(coms, disase_set):
    com_counts = [len(disase_set.intersection(set(coms[key]))) for key in coms.keys()]
    non_zero_com_counts = [count for count in com_counts if count > 0]
    non_zero_com_counts.sort(reverse=True)
    return non_zero_com_counts

def get_size_of_com_with_term(coms, disase_set, key_term):
    for key in coms.keys():
        if key_term in coms[key]:
            print(len(disase_set.intersection(set(coms[key]))))
            return disase_set.intersection(set(coms[key]))
    return None


print('AA')
print('Markov', get_non_zero_com_counts(m_coms, aa))
print('Louvain', get_non_zero_com_counts(l_coms, aa))
print('Walktrap', get_non_zero_com_counts(w_coms, aa))
print()
print('TOF')
print('Markov', get_non_zero_com_counts(m_coms, tof))
print('Louvain', get_non_zero_com_counts(l_coms, tof))
print('Walktrap', get_non_zero_com_counts(w_coms, tof))
print()
print('CFTR')
print('Markov', get_non_zero_com_counts(m_coms, cftr))
print('Louvain', get_non_zero_com_counts(l_coms, cftr))
print('Walktrap', get_non_zero_com_counts(w_coms, cftr))
print('Markov', get_size_of_com_with_term(m_coms, cftr,'CFTR'))
print('Louvain', get_size_of_com_with_term(l_coms, cftr,'CFTR'))
print('Walktrap', get_size_of_com_with_term(w_coms, cftr,'CFTR'))
print()
print('HTT')
print('Markov', get_non_zero_com_counts(m_coms, htt))
print('Louvain', get_non_zero_com_counts(l_coms, htt))
print('Walktrap', get_non_zero_com_counts(w_coms, htt))
print('Markov', get_size_of_com_with_term(m_coms, htt,'HTT'))
print('Louvain', get_size_of_com_with_term(l_coms, htt,'HTT'))
print('Walktrap', get_size_of_com_with_term(w_coms, htt,'HTT'))
