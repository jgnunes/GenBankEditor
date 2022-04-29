from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import sys

def get_feat_info(feat):

    feat_start = feat.location.start
    feat_end = feat.location.end
    feat_strand = feat.location.strand
    feat_len = len(feat.location) 
     
    if "product" in feat.qualifiers:
        feat_name = feat.qualifiers["product"][0]
    elif "gene" in feat.qualifiers:
        feat_name = feat.qualifiers["gene"][0]
    else:
        print(f"Error. Couldn't find feature name for feature {feat}.")
        exit(1)

    return (feat_start, feat_end, feat_strand, feat_len, feat_name) 

def rotate_genbank(in_gbk, ref_gene, out_gbk):
    
    for seq_record in SeqIO.parse(in_gbk, "genbank"):
        complete_seq_len = len(seq_record.seq)
        for feat in seq_record.features:
            if feat.type == "gene":
                feat_start, feat_end, feat_strand, feat_len, feat_name = get_feat_info(feat)
                if feat_name == ref_gene:
                    ref_feat = feat
                    ref_start, ref_end, ref_strand, ref_len, ref_name = feat_start, feat_end, feat_strand, feat_len, feat_name
    
    if ref_start and ref_end and ref_strand and ref_len:
        #print(ref_len)
        pass 
    else:
        print("Reference gene not found")
        sys.exit(1)

    with open(out_gbk, "w") as f: 
        for seq_record in SeqIO.parse(in_gbk, "genbank"):
            seq_record.id = f"{seq_record.id}_rotated"
            seq_record.name = f"{seq_record.name}_rotated"
            seq_record.description = ""
            # rotate feature coordinates
            for feat in seq_record.features:

                feat_start, feat_end, feat_strand, feat_len, feat_name = get_feat_info(feat)
                
                if feat_name == ref_name:
                    new_start = 0
                    new_end = new_start + feat_len
                    feat.location = FeatureLocation(start=new_start, end=new_end, strand=feat_strand)

                elif feat_start > ref_start:
                    new_start = feat_start - ref_start
                    new_end = new_start + feat_len
                    feat.location = FeatureLocation(start=new_start, end=new_end, strand=feat_strand)
                
                # else if feature comes befores the reference gene 
                elif feat_start < ref_start:
                    if feat_end < ref_start:
                        new_start = ref_len + (complete_seq_len - ref_end) + feat_start 
                        new_end = new_start + feat_len 
                        feat.location = FeatureLocation(start=new_start, end=new_end, strand=feat_strand)
                    else:
                        print("Warning: feature will possibly be truncated and divided into two features (before and after ORIGIN):")
                        print(feat)
                         

#                        # feature before ORIGIN
#                        new_start = ref_len + (complete_seq_len - ref_end) + feat_start
#                        new_end = new_start + feat_len
#                        feat.location = FeatureLocation(start=new_start, end=new_end, strand=feat_strand)
#                        if "product" in feat.qualifiers:
#                            feat.qualifiers["product"] = feat.qualifiers["product"][0] + "_beforeOrigin"
#                        elif "gene" in feat.qualifiers:
#                            feat.qualifiers["gene"] = feat.qualifiers["gene"][0] + "_beforeOrigin"
#                        
#                        # feature after ORIGIN
#                        overlap_len = feat_end - ref_start
#                        print(f"overlap_len: {overlap_len}")
#                        my_feat_location = FeatureLocation(start=0, end=overlap_len, strand=feat_strand)
#                        my_feat_type = feat.type
#                        
#                        if "product" in feat.qualifiers:
#                            #print("type of feat qualifiers product: ", type(feat.qualifiers["product"]))
#                            my_feat_name = feat.qualifiers["product"].replace("_beforeOrigin", "_afterOrigin")
#                            #print(f"my_feat_name: {my_feat_name}")
#                            #print(f"my_feat_location: {my_feat_location}")
#                            my_feature = SeqFeature(location=my_feat_location, type=my_feat_type, strand=feat.strand)
#                            my_feature.qualifiers["product"] = [my_feat_name]
#                            print("my_feature:")
#                            print(my_feature)
#                            seq_record.features.append(SeqFeature(location=my_feat_location, type=my_feat_type, strand=feat.strand))
#                        elif "gene" in feat.qualifiers:
#                            #my_feat_name = feat.qualifiers["gene"][0] + "_afterOrigin"
#                            my_feature = SeqFeature(location=my_feat_location,type=my_feat_type, strand=feat.strand)
#                            my_feature.qualifiers["gene"] = feat.qualifiers["gene"][0] + "_afterOrigin" 
#
#                        #seq_record.features.append(SeqFeature(my_feat_location)) 

            # rotate record sequence
            before_ref_start = seq_record.seq[0:ref_start]
            after_ref_start = seq_record.seq[ref_start:]
            rotated_seq = after_ref_start + before_ref_start
            seq_record.seq = rotated_seq
            
            # write rotated output files 
            SeqIO.write(seq_record, f, "genbank")
    
    # generate rotated fasta file
    SeqIO.convert(out_gbk, "genbank", f"{out_gbk}.fasta", "fasta")

if __name__ == "__main__":
        
    if sys.argv[1] == "-h":
        print("""python rotate_genbank.py <in_gbk> <ref_gene> <out_gbk>            

    in_gbk: input genbank file to be modified
    species_name: name of gene to be used as reference for rotation
    out_gbk: name of new, modified genbank file
                """)
        sys.exit(0)
    
    in_gbk = sys.argv[1]
    ref_gene = sys.argv[2]
    out_gbk = sys.argv[3]

    rotate_genbank(in_gbk, ref_gene, out_gbk)
