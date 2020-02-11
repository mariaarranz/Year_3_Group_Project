import csv
from decimal import Decimal

#------------------------------------------------------------------------------------------
#THIS CODE WILL OPEN THE FILE DATABASE AND CUT OUT THE SEQUENCES
with open('GS_mRNA.fsa') as database:
    trash = ''
    sequence = []
    i = 0
    for line in database.readlines():
        if line[0] == '>':
            trash += line
            if len(sequence) == 0:
                i = 0
            else:
                i += 1
        else:
            if len(sequence) == i:
                sequence.append(line)
            else:
                sequence[i] += line

#length_sequences = len(sequence)
#seq_num = input('Select a sequence from 0 to ' + str(length_sequences) + ':')
#n = int(seq_num)

n = 0
first_seq = sequence[n]
#print(first_seq)



#------------------------------------------------------------------------------------------
#THIS WILL TRANSLATE THE DNA SEQUENCE PICKED (HERE SEQUENCE 0) INTO RNA SEQUENCE
trna_seq = ''
for i in first_seq:
    if i == 'A':
        trna_seq += 'U'
    elif i == 'T':
        trna_seq += 'A'
    elif i == 'C':
        trna_seq += 'G'
    elif i == 'G':
        trna_seq += 'C'

#print(trna_seq)


#------------------------------------------------------------------------------------------
#THIS WILL CREATE A LIST OF ALL CODONS BY SEPARATING THE SEQUENCE INTO GROUPS OF 3
codons = []
start = 0

for l in range(int(len(trna_seq)/3)):
    end = start + 3
    codon = trna_seq[start:end]
    codons.append(codon)
    start = end

#print(codons)


#------------------------------------------------------------------------------------------
#RELATIVE SYNONYMOUS CODON USAGE CALCULATION
translation_list = {'*':['UAA','UAG','UGA'],
                    'A':['GCA','GCC','GCG','GCU'],
                    'C':['UGC','UGU'],
                    'D':['GAC','GAU'],
                    'E':['GAA','GAG'],
                    'F':['UUC','UUU'],
                    'G':['GGA','GGC','GGG','GGU'],
                    'H':['CAC','CAU'],
                    'I':['AUA','AUC','AUU'],
                    'K':['AAA','AAG'],
                    'L':['CUA','CUC','CUG','CUU','UUA','UUG'],
                    'M':['AUG'],
                    'N':['AAC','AAU'],
                    'P':['CCA','CCC','CCG','CCU'],
                    'Q':['CAA','CAG'],
                    'R':['AGA','AGG','CGA','CGC','CGG','CGU'],
                    'S':['AGC','AGU','UCA','UCC','UCG','UCU'],
                    'T':['ACA','ACC','ACG','ACU'],
                    'V':['GUA','GUC','GUG','GUU'],
                    'W':['UGG'],
                    'Y':['UAC','UAU']}

#RSCU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
count_dic = []
rscu = []
for aa in translation_list:
    total = 0.0
    s_codons = translation_list[aa]

    #THIS WILL COUNT HOW MANY TIMES EACH CODON APPEARS IN THE SEQUENCE
    for s in s_codons:
        count = codons.count(s)
        count_dic.append(count)
        total += count

    #calculate the RSCU value for each of the codons
    for  s in s_codons:
        count = codons.count(s)
        denominator = float(total)/ len(s_codons)
        r = Decimal(count/denominator)
        result = round(r,2)
        rscu.append(result)     


#------------------------------------------------------------------------------------------
#THIS WILL TRANSLATE THE CODONS INTO THE CORRESPONDING AMINO ACIDS
translation = ''
for c in codons:
    for aa in translation_list:
        s_codons = translation_list[aa]
        for s in s_codons:
            if s == c:
                translation += aa

#print(translation)           


#------------------------------------------------------------------------------------------
#THIS WILL SPLIT THE TRANSLATED AMINO ACIDS BY ITS STOP CODONS

split_stop = translation.split('*')
#print(split_stop)


#THIS WILL CREATE A LIST WITH THE SECTIONS THAT CONTAIN THE START CODON
s = 'M'
proteins = []
for a in split_stop:
    if s in a:
        proteins.append(a)

#print(proteins)


#THIS WILL THEN FINISH SPLITTING THE PROTEINS FROM START CODON TO STOP CODON
final_p = []
for p in proteins:
    start = p.find(s)
    complete = p[start:]
    final_p.append(complete)

#print(final_p)
    

#------------------------------------------------------------------------------------------
#THIS WILL CREATE AN OUTPUT EXCEL FILE WITH THE FREQUENCY TABLE
excel_list = [{'Amino Acid':'*', 'Codon':'UAA', 'Count':count_dic[0], 'RSCU': rscu[0]},
              {'Amino Acid':'*', 'Codon':'UAG', 'Count':count_dic[1], 'RSCU': rscu[1]},
              {'Amino Acid':'*', 'Codon':'UGA', 'Count':count_dic[2], 'RSCU': rscu[2]},
              {'Amino Acid':'A', 'Codon':'GCA', 'Count':count_dic[3], 'RSCU': rscu[3]},
              {'Amino Acid':'A', 'Codon':'GCC', 'Count':count_dic[4], 'RSCU': rscu[4]},
              {'Amino Acid':'A', 'Codon':'GCG', 'Count':count_dic[5], 'RSCU': rscu[5]},
              {'Amino Acid':'A', 'Codon':'GCU', 'Count':count_dic[6], 'RSCU': rscu[6]},
              {'Amino Acid':'C', 'Codon':'UGC', 'Count':count_dic[7], 'RSCU': rscu[7]},
              {'Amino Acid':'C', 'Codon':'UGU', 'Count':count_dic[8], 'RSCU': rscu[8]},
              {'Amino Acid':'D', 'Codon':'GAC', 'Count':count_dic[9], 'RSCU': rscu[9]},
              {'Amino Acid':'D', 'Codon':'GAU', 'Count':count_dic[10], 'RSCU': rscu[10]},
              {'Amino Acid':'E', 'Codon':'GAA', 'Count':count_dic[11], 'RSCU': rscu[11]},
              {'Amino Acid':'E', 'Codon':'GAG', 'Count':count_dic[12], 'RSCU': rscu[12]},
              {'Amino Acid':'F', 'Codon':'UUC', 'Count':count_dic[13], 'RSCU': rscu[13]},
              {'Amino Acid':'F', 'Codon':'UUU', 'Count':count_dic[14], 'RSCU': rscu[14]},
              {'Amino Acid':'G', 'Codon':'GGA', 'Count':count_dic[15], 'RSCU': rscu[15]},
              {'Amino Acid':'G', 'Codon':'GGC', 'Count':count_dic[16], 'RSCU': rscu[16]},
              {'Amino Acid':'G', 'Codon':'GGG', 'Count':count_dic[17], 'RSCU': rscu[17]},
              {'Amino Acid':'G', 'Codon':'GGU', 'Count':count_dic[18], 'RSCU': rscu[18]},
              {'Amino Acid':'H', 'Codon':'CAC', 'Count':count_dic[19], 'RSCU': rscu[19]},
              {'Amino Acid':'H', 'Codon':'CAU', 'Count':count_dic[20], 'RSCU': rscu[20]},
              {'Amino Acid':'I', 'Codon':'AUA', 'Count':count_dic[21], 'RSCU': rscu[21]},
              {'Amino Acid':'I', 'Codon':'AUC', 'Count':count_dic[22], 'RSCU': rscu[22]},
              {'Amino Acid':'I', 'Codon':'AUU', 'Count':count_dic[23], 'RSCU': rscu[23]},
              {'Amino Acid':'K', 'Codon':'AAA', 'Count':count_dic[24], 'RSCU': rscu[24]},
              {'Amino Acid':'K', 'Codon':'AAG', 'Count':count_dic[25], 'RSCU': rscu[25]},
              {'Amino Acid':'L', 'Codon':'CUA', 'Count':count_dic[26], 'RSCU': rscu[26]},
              {'Amino Acid':'L', 'Codon':'CUC', 'Count':count_dic[27], 'RSCU': rscu[27]},
              {'Amino Acid':'L', 'Codon':'CUG', 'Count':count_dic[28], 'RSCU': rscu[28]},
              {'Amino Acid':'L', 'Codon':'CUU', 'Count':count_dic[29], 'RSCU': rscu[29]},
              {'Amino Acid':'L', 'Codon':'UUA', 'Count':count_dic[30], 'RSCU': rscu[30]},
              {'Amino Acid':'L', 'Codon':'UUG', 'Count':count_dic[31], 'RSCU': rscu[31]},
              {'Amino Acid':'M', 'Codon':'AUG', 'Count':count_dic[32], 'RSCU': rscu[32]},
              {'Amino Acid':'N', 'Codon':'AAC', 'Count':count_dic[33], 'RSCU': rscu[33]},
              {'Amino Acid':'N', 'Codon':'AAU', 'Count':count_dic[34], 'RSCU': rscu[34]},
              {'Amino Acid':'P', 'Codon':'CCA', 'Count':count_dic[35], 'RSCU': rscu[35]},
              {'Amino Acid':'P', 'Codon':'CCC', 'Count':count_dic[36], 'RSCU': rscu[36]},
              {'Amino Acid':'P', 'Codon':'CCG', 'Count':count_dic[37], 'RSCU': rscu[37]},
              {'Amino Acid':'P', 'Codon':'CCU', 'Count':count_dic[38], 'RSCU': rscu[38]},
              {'Amino Acid':'Q', 'Codon':'CAA', 'Count':count_dic[39], 'RSCU': rscu[39]},
              {'Amino Acid':'Q', 'Codon':'CAG', 'Count':count_dic[40], 'RSCU': rscu[40]},
              {'Amino Acid':'R', 'Codon':'AGA', 'Count':count_dic[41], 'RSCU': rscu[41]},
              {'Amino Acid':'R', 'Codon':'AGG', 'Count':count_dic[42], 'RSCU': rscu[42]},
              {'Amino Acid':'R', 'Codon':'CGA', 'Count':count_dic[43], 'RSCU': rscu[43]},
              {'Amino Acid':'R', 'Codon':'CGC', 'Count':count_dic[44], 'RSCU': rscu[44]},
              {'Amino Acid':'R', 'Codon':'CGG', 'Count':count_dic[45], 'RSCU': rscu[45]},
              {'Amino Acid':'R', 'Codon':'CGU', 'Count':count_dic[46], 'RSCU': rscu[46]},
              {'Amino Acid':'S', 'Codon':'AGC', 'Count':count_dic[47], 'RSCU': rscu[47]},
              {'Amino Acid':'S', 'Codon':'AGU', 'Count':count_dic[48], 'RSCU': rscu[48]},
              {'Amino Acid':'S', 'Codon':'UCA', 'Count':count_dic[49], 'RSCU': rscu[49]},
              {'Amino Acid':'S', 'Codon':'UCC', 'Count':count_dic[50], 'RSCU': rscu[50]},
              {'Amino Acid':'S', 'Codon':'UCG', 'Count':count_dic[51], 'RSCU': rscu[51]},
              {'Amino Acid':'S', 'Codon':'UCU', 'Count':count_dic[52], 'RSCU': rscu[52]},
              {'Amino Acid':'T', 'Codon':'ACA', 'Count':count_dic[53], 'RSCU': rscu[53]},
              {'Amino Acid':'T', 'Codon':'ACC', 'Count':count_dic[54], 'RSCU': rscu[54]},
              {'Amino Acid':'T', 'Codon':'ACG', 'Count':count_dic[55], 'RSCU': rscu[55]},
              {'Amino Acid':'T', 'Codon':'ACU', 'Count':count_dic[56], 'RSCU': rscu[56]},
              {'Amino Acid':'V', 'Codon':'GUA', 'Count':count_dic[57], 'RSCU': rscu[57]},
              {'Amino Acid':'V', 'Codon':'GUC', 'Count':count_dic[58], 'RSCU': rscu[58]},
              {'Amino Acid':'V', 'Codon':'GUG', 'Count':count_dic[59], 'RSCU': rscu[59]},
              {'Amino Acid':'V', 'Codon':'GUU', 'Count':count_dic[60], 'RSCU': rscu[60]},
              {'Amino Acid':'W', 'Codon':'UGG', 'Count':count_dic[61], 'RSCU': rscu[61]},
              {'Amino Acid':'Y', 'Codon':'UAC', 'Count':count_dic[62], 'RSCU': rscu[62]},
              {'Amino Acid':'Y', 'Codon':'UAU', 'Count':count_dic[63], 'RSCU': rscu[63]}]



amino_dict = {'A':'alanine',                'P':'proline',
              'B':'aspartate/asparagine',   'Q':'glutamine',
              'C':'cystine',                'R':'arginine',
              'D':'aspartate',              'S':'serine',
              'E':'glutamate',              'T':'threonine',
              'F':'phenylalanine',          'U':'selenocysteine',
              'G':'glycine',                'V':'valine',
              'H':'histidine',              'W':'tryptophan',
              'I':'isoleucine',             'Y':'tyrosine',
              'K':'lysine',                 'Z':'glutamate/glutamine',
              'L':'leucine',                '*':'STOP',
              'M':'methionine',             
              'N':'asparagine'}              


with open('output.csv', 'w') as output_csv:
    fields = ['Amino Acid', 'Codon', 'Count', 'RSCU']
    output_writer = csv.DictWriter(output_csv, fieldnames=fields)

    output_writer.writeheader()
    for item in excel_list:
        output_writer.writerow(item)
