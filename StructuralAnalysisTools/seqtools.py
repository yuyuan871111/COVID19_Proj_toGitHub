import re

def seq2seqList(seq):
    seqList = []
    for each_aa in seq:
        seqList.append(each_aa)
    return seqList

def seqList2seq(seqList):
    return ''.join(seqList)

def check_mismatch(test1,test2, calib_number = 0):
    mismatch = False
    for index, (i, j) in enumerate(zip(test1, test2)):
        if not i == j:
            print(f'{index+1+calib_number}: mismatch -- {test1[index]} to {test2[index]}')
            mismatch = True
        if index+1 == len(test1) and mismatch == False:
            print("No mismatch.")

def create_variants_seq(variants, ref_seq):
    # input variants => e.g. 'T19R, (G142D), 156del, 157del, R158G, L452R, T478K, D614G, P681R, D950N'
    variants = variants.split(sep=', ')
    variants_points = []
    for variants_point in variants:
        m = re.search('(.\d+)(\w+)', variants_point)
        if m.group(2) == 'del':
            variants_points.append([m.group(2), m.group(1), ''])
        else:
            m_ = re.search('(\w)(\d+)', m.group(1))
            variants_points.append([m_.group(1), m_.group(2), m.group(2)])   
    variants_seq = ref_seq.copy()
    for i, variants_point in enumerate(variants_points):
        pos = int(variants_point[1])
        if not variants_point[0] == 'del':
            check_mismatch(variants_point, variants_seq[pos-1])
            print(f'variants checked...{i+1}')
        variants_seq[pos-1] = variants_point[2]
    return variants_seq

def read_fasta_seq(file_path):
    with open(file_path, 'r') as F:
        fasta_data = F.read().splitlines()
    header_index = []
    for index, eachline in enumerate(fasta_data):
        m = re.search('^>', eachline)
        if not m == None:
            header_index.append(index)
    
    if len(header_index) == 1:
        return seq2seqList(seqList2seq(fasta_data[1:]))
    else:
        data_dict = {}
        for idx_in_list, index in enumerate(header_index):
            if not (idx_in_list+1) == len(header_index):
                data_dict[fasta_data[index]] = seq2seqList(seqList2seq(fasta_data[index+1 : header_index[idx_in_list+1]]))
            else:
                data_dict[fasta_data[index]] = seq2seqList(seqList2seq(fasta_data[index+1:]))
        return data_dict
