from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

if sys.argv[1] == "-h":
    print("""python add_organism.py <in_gbk> <species_name> <out_gbk>(optional)            

in_gbk: input genbank file to be modified
species_name: name of species to be attributed
out_gbk: name of new, modified genbank file (default="new_gb.bg")
            """)
    sys.exit(0)

in_gbk = sys.argv[1]
organism = sys.argv[2]

if len(sys.argv)==4:
    out_file = sys.argv[3]
else:
    out_file = "new_gb.gb"

with open(out_file, "w") as f:
    for seq_record in SeqIO.parse(in_gbk, "genbank"):
        seq_len = len(seq_record.seq)
        # add organism annotation
        seq_record.annotations['organism'] = organism
        # add feature with organism qualifier
        # create source feature
        my_feat_location = FeatureLocation(0, seq_len)
        my_feat_type = "source"
        my_feature = SeqFeature(my_feat_location, type=my_feat_type)
        my_feature.qualifiers['organism'] = organism
        seq_record.features.append(my_feature)
        SeqIO.write(seq_record, f, "genbank")
