import hashlib

def convert_psl_ids(psl_file, mapping_file):
    """
    Convert read names in a PSL file to stable hash-based IDs 
    based on splice junction starts and ends.

    Parameters:
        psl_file (str): Input PSL file
        mapping_file (str): Output file for ID mappings
    """
    mappings = []
    
    with open(psl_file, "r") as infile:
        for line in infile:
            read          = line.rstrip("\n").split("\t")
            original_name = read[9]
            chrom         = read[13]
            strand        = read[8]

            # parse blocksizes, tstarts
            blocksizes = [int(x) for x in read[18].rstrip(",").split(",")]
            tstarts    = [int(x) for x in read[20].rstrip(",").split(",")]
            
            # Check for monoexonic transcripts (single exon = no junctions)
            if len(blocksizes) == 1:
                new_id = f"{chrom}_monoexon:{tstarts[0]}-{tstarts[0] + blocksizes[0]}_{strand}"
            else:
                # compute tends
                tends = [s + l for s, l in zip(tstarts, blocksizes)]

                # collapse junctions
                collapse_tstarts = tstarts[1:]
                collapse_tends   = tends[:-1]

                # hash start junctions
                s        = ",".join(map(str, collapse_tstarts))
                s_hashed = hashlib.shake_256(s.encode("utf-8")).hexdigest(8)
                s_id     = "s" + s_hashed

                # hash end junctions
                e        = ",".join(map(str, collapse_tends))
                e_hashed = hashlib.shake_256(e.encode("utf-8")).hexdigest(8)
                e_id     = "e" + e_hashed

                # new read ID = chr_startID:endID_strand
                new_id  = f"{chrom}_{s_id}:{e_id}_{strand}"
                
            mappings.append((original_name, new_id))
            
    # Write mapping file
    with open(mapping_file, "w") as mapfile:
        mapfile.write("transcript_id\tgene_id\thash_id\n")
        for original, hashed in mappings:
          
            # Split transcript_id_gene_id format
            parts = original.split("_")
            if len(parts) == 2:
                transcript_id, gene_id = parts
                mapfile.write(f"{transcript_id}\t{gene_id}\t{hashed}\n")


# Run function
if __name__ == "__main__":
    import sys
    psl_in = sys.argv[1]
    mapping_out = sys.argv[2]
    convert_psl_ids(psl_in, mapping_out)
