import pysam
import numpy as np
import matplotlib.pyplot as plt
import glob, os

plt.close('all')

#%%

path = os.getcwd()

tp_bwa, fp_bwa, fn_bwa = [], [], []
tp_bow, fp_bow, fn_bow = [], [], []

for file in glob.glob(os.path.join(path, 'data1_*.txt')):
    
    sett = file[-5]
    f = open(file, "r")
    lines = f.readlines()
    
    bwa_pct = lines[1]
    bwa_pct = bwa_pct.split()
    bwa_pct = [float(i) for i in bwa_pct]
    tp_bwa.append(bwa_pct[0])
    fp_bwa.append(bwa_pct[1])
    fn_bwa.append(bwa_pct[2])
    
    bow_pct = lines[9]
    bow_pct = bow_pct.split()
    bow_pct = [float(i) for i in bow_pct]
    tp_bow.append(bow_pct[0])
    fp_bow.append(bow_pct[1])
    fn_bow.append(bow_pct[2])
    
    bwa_mapq_a = lines[4]
    bwa_mapq_a = bwa_mapq_a.split()
    bwa_mapq_a = [float(i) for i in bwa_mapq_a]
    bwa_mapq_m = lines[6]
    bwa_mapq_m = bwa_mapq_m.split()
    bwa_mapq_m = [float(i) for i in bwa_mapq_m]
    
    bow_mapq_a = lines[12]
    bow_mapq_a = bow_mapq_a.split()
    bow_mapq_a = [float(i) for i in bow_mapq_a]
    bow_mapq_m = lines[14]
    bow_mapq_m = bow_mapq_m.split()
    bow_mapq_m = [float(i) for i in bow_mapq_m]
    
    data = [bwa_mapq_a, bwa_mapq_m, bow_mapq_a, bow_mapq_m]
    sets = np.linspace(1, 6, 6)
    label = ['bwa-mem aligned', 'bwa-mem misaligned', 'bowtie aligned', \
             'bowtie misaligned']
    fig = plt.figure(figsize=(12,8))
    for c, vals in enumerate(data):
        if c == 0 or c == 2:
            pp = plt.bar(sets, vals, width=0.45, label=label[c])
        if c == 1:
            pp = plt.bar(sets, vals, width=0.45, bottom=data[0], \
                    label=label[c])
        if c == 3:
            pp = plt.bar(sets, vals, width=0.45, bottom=data[2], \
                    label=label[c])
        plt.bar_label(pp, label_type='center', size=10)
    plt.xticks([r+1 for r in range(6)],  ['0-10', '10-20', '20-30', \
                                                '30-40', '40-50', '50-60'])
    plt.title('Mapping quality distribution - set ' + sett)
    plt.ylabel('Percentage of reads')
    # plt.ylim((1,110))
    plt.legend()
    f.close()

data = [tp_bwa, tp_bow, fp_bwa, fp_bow, fn_bwa, fn_bow]
label = ['bwa-mem aligned', 'bowtie aligned', 'bwa-mem misaligned', \
          'bowtie misaligned', 'bwa-mem unmapped', 'bowtie unmapped']
sets = np.linspace(1, 3, 3)

fig = plt.figure(figsize=(14, 6))
for c, vals in enumerate(data):
    pp = plt.bar(sets + c*0.15, vals, width=0.15, label=label[c])
    plt.bar_label(pp, label_type='edge', size=8)

plt.xticks([r + 0.25 for r in [1, 2, 3]],  ['Set 1', 'Set 2', 'Set 3'])
plt.ylabel('Percentage of reads')
plt.legend(loc='best')

np_data = []
for arr in data:
    np_data.append(np.array(arr))
    
precision_bwa = np_data[0] / (np_data[0] + np_data[2])
precision_bow = np_data[1] / (np_data[1] + np_data[3])
fig = plt.figure(figsize=(14, 8))
pp = plt.bar(sets, precision_bwa, width=0.15, label='bwa')
plt.bar_label(pp, label_type='edge', size=8)
pp = plt.bar(sets + 0.15, precision_bow, width=0.15, label='bow')
plt.bar_label(pp, label_type='edge', size=8)
plt.xticks([r + 0.15/2 for r in [1, 2, 3]],  ['Set 1', 'Set 2', 'Set 3'])
plt.ylabel('Precision')
plt.legend(loc='best')

recall_bwa = np_data[0] / (np_data[0] + np_data[4])
recall_bow = np_data[1] / (np_data[1] + np_data[5])
fig = plt.figure(figsize=(14, 8))
pp = plt.bar(sets, recall_bwa, width=0.15, label='bwa')
plt.bar_label(pp, label_type='center', size=8)
pp = plt.bar(sets + 0.15, recall_bow, width=0.15, label='bow')
plt.bar_label(pp, label_type='center', size=8)
plt.xticks([r + 0.15/2 for r in [1, 2, 3]],  ['Set 1', 'Set 2', 'Set 3'])
plt.ylabel('Recall')
plt.legend(loc='best')

f_score_bwa = 2 * precision_bwa * recall_bwa / (precision_bwa + recall_bwa)
f_score_bow = 2 * precision_bow * recall_bow / (precision_bow + recall_bow)
fig = plt.figure(figsize=(14, 8))
pp = plt.bar(sets, f_score_bwa, width=0.15, label='bwa')
plt.bar_label(pp, label_type='center', size=8)
pp = plt.bar(sets + 0.15, f_score_bow, width=0.15, label='bow')
plt.bar_label(pp, label_type='center', size=8)
plt.xticks([r + 0.15/2 for r in [1, 2, 3]],  ['Set 1', 'Set 2', 'Set 3'])
plt.ylabel('F-score')
plt.legend(loc='best')