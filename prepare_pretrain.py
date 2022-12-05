import re
import os, glob
import pysam
import h5py
from tqdm import tqdm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

bin_size = 128; bin_num = 4; seq_overlap = 0
seq_len = bin_size * bin_num; shift = int(seq_len * (1 - seq_overlap))
bin_size, seq_len, shift

os.chdir("/public/home/guogjgroup/ggj/JiaqiLi/Nvwa2/0_Dataset")

genome = pysam.FastaFile("../../Resource/STAR_Reference_Mouse/Mus_musculus.GRCm38.88.fasta")
genome

references = genome.references[:22]
print(references)

seq_info = []

for chrom in tqdm(references):
    chrom_size = genome.get_reference_length(chrom)
    for start in range(0, chrom_size, shift):
        end = start + seq_len
        if end > chrom_size:
            pad = 'N' * (end - chrom_size) # pad N
            end = chrom_size
            # print("padding chrom%s:%d-%d"%(chrom, start, end))
    
        info = "%s:%d-%d" % (chrom, start, end)
        seq_info.append([info, chrom, start, end])
    
seq_info = np.asarray(seq_info)
genome.close()

seq_info.shape

seq_info = np.hstack([np.zeros((seq_info.shape[0], 1)), seq_info])
seq_info[:,0] = range(seq_info.shape[0])

files = sorted(glob.glob("../../Resource/Encode_mm10/encode_mm10/*.sort.bed.gz"))[663:]
len(files)

def create_npz(fname):
    tabixfile = pysam.TabixFile(fname)
    labels_anno = fname.split("/")[-1]

    label = []
    chroms_except = []
    for chrom in references:
        chrom_size = genome.get_reference_length(chrom)
        for start in range(0, chrom_size, shift):
            end = start + seq_len
            if end > chrom_size:
                # print("padding chrom%s:%d-%d"%(chrom, start, end))
                end = chrom_size

            peak_label = np.zeros(bin_num, dtype=np.float32)
            try:
                for line in tabixfile.fetch("chr"+chrom, start, end, parser=pysam.asTuple()):
                    _, peak_start, peak_end = line[:3]
                    peak_start, peak_end = int(peak_start) - start, int(peak_end) - start # relative position
                    peak_start, peak_end = peak_start // bin_size, peak_end // bin_size + 1
                    peak_label[peak_start:peak_end] = 1. # replace with positive labels
            except:
                chroms_except.append(chrom)
            label.append(peak_label)

    tabixfile.close()
    label = np.vstack(label)[:, :, None] # unsqueeze_aixs=-1 
    chroms_except = ";".join(np.unique(chroms_except).astype(str))

    fname = os.path.join("../../Resource/Encode_mm10/encode_mm10_npz", labels_anno.replace("gz", "npz"))
    np.savez_compressed(fname, label=label, chroms_except=chroms_except)

from multiprocessing import Pool
pool = Pool()
pool.map(create_npz, files)
pool.close()
pool.join()

# for f in files:
#     create_npz(f)


files = sorted(glob.glob("../../Resource/Encode_mm10/encode_mm10_npz/*.sort.bed.npz"))
len(files)

labels, labels_anno, chroms_except = [], [], []

for fname in tqdm(files):
    npz = np.load(fname)
    species, dtype = fname.split("/")[-2:]
    dtype = dtype.split(".")[0]
    labels_anno.append([species, dtype])
    labels.append(npz['label'])
    chroms_except.append(npz['chroms_except'])

labels = np.concatenate(labels, axis=-1)
labels_anno = np.array(labels_anno)
labels.shape, labels_anno.shape

label_regulatory = labels.copy()
label_regulatory_anno = labels_anno.copy()



output_dir = "Dataset.mm10.Encode.h5"
data_type = np.float32 # bool
compress_args = {'compression': 'gzip', 'compression_opts': 1}

h5file = h5py.File(output_dir, 'r+')
print(h5file.keys())
del h5file["label_regulatory"]
del h5file["label_regulatory_anno"]
h5file.create_dataset("label_regulatory", data=label_regulatory, dtype=data_type, **compress_args)
h5file.create_dataset("label_regulatory_anno", data=label_regulatory_anno.astype(np.bytes_), **compress_args)
h5file.close()

h5file = h5py.File(output_dir, 'r')
print(h5file.keys())
h5file.close()