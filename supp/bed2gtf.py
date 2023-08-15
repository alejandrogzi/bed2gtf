#!/usr/bin/env python3



import time
import argparse
import sys



__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"



GENES = {}


def check_exon_frame(exonStart, exonEnd, cdsStart, cdsEnd):
    start = max(exonStart, cdsStart)
    end = min(exonEnd, cdsEnd)
    return start, end



def hashmap(path):
    file = open(path, "r")
    isoforms = {}
    for line in file:
        if line.startswith('Gene'):
            continue
        else:
            fields = line.strip().split('\t')
            isoforms[fields[1].rstrip()] = fields[0]
    return isoforms



def exon_frames(cdsStart: int, cdsEnd: int, exonStart: list, exonEnd: list, strand: str):
    frames = [0]*len(exonStart)
    cds = 0
    start = 0 if strand == "+" else len(exonStart) - 1
    end = len(exonStart) if strand == "+" else -1 
    step = 1 if strand == "+" else -1
    
    for exon in range(start, end, step):
        cdsExonStart, cdsExonEnd = check_exon_frame(exonStart[exon], exonEnd[exon], cdsStart, cdsEnd)
        if cdsExonStart < cdsExonEnd:
            frames[exon] = cds % 3
            cds += (cdsExonEnd - cdsExonStart)
        else:
            frames[exon] = -1
            
    return frames 



def find_first_codon(frames, strand, cdsStart, cdsEnd, exonStart, exonEnd):
    codon = [0]*6
    exon_index = 0
    for fr in range(len(frames)-1): #extracts the index of the first exon that has a valid CDS 
        if frames[fr] >= 0:
            exon_index = fr
            break
        else:
            return codon
    
    cdsExonStart, cdsExonEnd = check_exon_frame(exonStart[exon_index], exonEnd[exon_index], cdsStart, cdsEnd)
    frame = frames[exon_index] if strand == "+" else (frames[exon_index] + (cdsExonEnd - cdsExonStart)) % 3

    if frame != 0: # not on a frame boundary
        return [0]*6
        
    codon[:2] = [cdsStart, cdsStart + min(cdsEnd - cdsStart, 3)]
    codon[4] = exon_index
    
    if codon[1] - codon[0] < 3: # checks if the start_codon is valid (must be 3 bases long)
            codon[5] = exon_index
            codon[4] = exon_index + 1
            codon[2] = cdsStart #the start will be always the same
            needed = 3 - (codon[1] - codon[0]) # calculates how many bases are left to complete a codon with 3 bases

            if cdsEnd - cdsStart < needed: # if the CDS is not long enough, then there is no valid start_codon (can't complete 3 bases long)
                return [0]*6

            codon[3] = codon[2] + needed # completes the bases left to achieve a valid start_codon
    
    return codon



def find_last_codon(frames, strand, cdsStart, cdsEnd, exonStart, exonEnd):
    codon = [0]*6
    exon_index = 0
    for fr in range(len(frames)-1,-1,-1):
        if frames[fr] >= 0:
            exon_index = fr
            break
        else:
            return [0]*6
    
    cdsExonStart, cdsExonEnd = check_exon_frame(exonStart[exon_index], exonEnd[exon_index], cdsStart, cdsEnd)
    frame = frames[exon_index] if strand == "-" else (frames[exon_index] + (cdsExonEnd - cdsExonStart)) % 3

    if frame != 0: # not on a frame boundary
        return codon
    
    codon[:2] = [max(cdsStart, cdsEnd-3), cdsEnd]
    codon[4] = exon_index

    if codon[1] - codon[0] < 3: # checks if the start_codon is valid (must be 3 bases long)
            codon[5] = exon_index
            codon[4] = exon_index - 1
            codon[2] = cdsStart #the start will be always the same
            needed = 3 - (codon[1] - codon[0])  # calculates how many bases are left to complete a codon with 3 bases

            if cdsEnd - cdsStart < needed: # if the CDS is not long enough, then there is no valid start_codon (can't complete 3 bases long)
                return [0]*6

            codon[3] = codon[2] + needed # completes the bases left to achieve a valid start_codon
    
    return codon



def codon_complete(codon):
    return (((codon[1] - codon[0]) + (codon[3]- codon[2])) == 3)



def in_exon(exon, pos, exonStarts, exonEnds):
    return exonStarts[exon] <= pos <= exonEnds[exon]



def move_pos(txStart, txEnd, exonCount, exonStart, exonEnd, pos, dist):
    orig_pos = pos
    assert int(txStart) <= pos <= int(txEnd)
    
    exon = None
    for i in range(exonCount):
        if in_exon(i, pos, exonStart, exonEnd): # checks if pos (cdsStart or cdsEnd) are within exon boundaries
            exon = i
            break
    
    if exon is None: # if pos is not within any exon boundary, return error
        err_msg = "can't find {}"
        err_msg = err_msg.format(pos)
        raise ValueError(err_msg)
    
    steps = abs(dist)
    direction = 1 if dist >= 0 else -1 # negative (backwards) for positive strand and viceversa
    while 0 <= exon < int(exonCount) and steps > 0:
        if in_exon(exon, pos + direction, exonStart, exonEnd): # if we sum -1/1 to the position (cdsStart or cdsEnd), are we still within exon boundaries? 
            pos += direction
            steps -= 1
        elif direction >= 0: # are we a negative strand?
            exon += 1
            if exon < exonCount:
                pos = exonStart[exon]
        else:
            exon -= 1
            if exon >= 0:
                pos = exonEnd[exon] - 1
                steps -= 1
    
    if steps > 0:
        err_msg = "can't move {} by {}"
        err_msg = err_msg.format(orig_pos, dist)
        raise ValueError(err_msg)
    
    return pos


def build_gene_line(source, geneName, chrom, strand, type, start, end, file):
    assert start <= end

    gtf_line = f"{chrom}\t"
    gtf_line += f"{source}\t"
    gtf_line += f"{type}\t"
    gtf_line += f"{start + 1}\t"
    gtf_line += f"{end}\t"
    gtf_line += ".\t" 
    gtf_line += f"{strand}\t"
    gtf_line += f".\t"
    gene_id = geneName
    gtf_line += f"gene_id \"{gene_id}\"; "
    gtf_line += "\n"
    file.write(gtf_line)

    return 



def build_gtf_line(source, name, chrom, strand, geneName, type_g, start, end, exon, frame, exonCount, file):
    assert start <= end

    # Convert frame to phase
    phase = '.' if frame < 0 else '0' if frame == 0 else '2' if frame == 1 else '1'
    gtf_line = f"{chrom}\t"
    gtf_line += f"{source}\t"
    gtf_line += f"{type_g}\t"
    gtf_line += f"{start + 1}\t"
    gtf_line += f"{end}\t"
    gtf_line += ".\t" 
    gtf_line += f"{strand}\t"
    gtf_line += f"{phase}\t"
    gene_id = geneName
    gtf_line += f"gene_id \"{gene_id}\"; "
    gtf_line += f"transcript_id \"{name}\"; "
    if exon >= 0:
        if strand == '-':
            gtf_line += f"exon_number \"{exonCount - exon}\"; "
            gtf_line += f"exon_id \"{name}.{exonCount - exon}\";"
        else:
            gtf_line += f"exon_number \"{exon + 1}\"; "
            gtf_line += f"exon_id \"{name}.{exon + 1}\";"
    gtf_line += "\n"
    file.write(gtf_line)

    return



def write_features(i, source, name, chrom, strand, geneName, firstUtrEnd, cdsStart, cdsEnd, lastUtrStart, frame, exonStart, exonEnd, exonCount, file):
    exonStart = exonStart[i]
    exonEnd = exonEnd[i]

    if exonStart < firstUtrEnd:
        end = min(exonEnd, firstUtrEnd)
        build_gtf_line(source, name, chrom, strand, geneName, "5UTR" if strand == "+" else "3UTR", exonStart, end, i, -1, exonCount, file)

    if cdsStart < exonEnd and cdsEnd > exonStart:
        start = max(exonStart, cdsStart)
        end = min(exonEnd, cdsEnd)
        build_gtf_line(source, name, chrom, strand, geneName, "CDS", start, end, i, frame, exonCount, file)

    if exonEnd > lastUtrStart:
        start = max(lastUtrStart, exonStart)
        build_gtf_line(source, name, chrom, strand, geneName, "3UTR" if strand == "+" else "5UTR", start, exonEnd, i, -1, exonCount, file)

    return



def write_codon(source, name, chrom, strand, geneName, type_g, codon, exonCount, file):
    build_gtf_line(source, name, chrom, strand, geneName, type_g, codon[0],codon[1], codon[4], 0, exonCount, file)
    
    if codon[2] < codon[3]:
        build_gtf_line(source, name, chrom, strand, geneName, type_g, codon[2],codon[1], codon[5], (codon[1] - codon[0]), exonCount, file)
    
    return



def run(line, file, isoforms):
    fields = line.strip().split('\t')
    #name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStart, exonEnd = fields #genePred
    chrom, txStart, txEnd, name, _, strand, cdsStart, cdsEnd, _, exonCount, exonEnd, exonStart = fields #bed

    txStart = int(txStart)
    txEnd = int(txEnd)
    cdsStart = int(cdsStart)
    cdsEnd = int(cdsEnd)
    exonCount = int(exonCount)
    geneName = isoforms[name]
    #exonStart = list(map(int, exonStart.split(',')[:-1])) #genePred
    #exonEnd = list(map(int, exonEnd.split(',')[:-1])) #genePred

    if exonCount != 0:
        exonStart = [txStart + x for x in list(map(int, exonStart.split(',')[:-1]))]
        exonEnd = [exonStart[i] + x for i,x in enumerate(list(map(int, exonEnd.split(',')[:-1])))]
    source = "test"
    
    frames = exon_frames(cdsStart, cdsEnd, exonStart, exonEnd, strand)
    
    first_codon = find_first_codon(frames, strand, cdsStart, cdsEnd, exonStart, exonEnd)
    last_codon = find_last_codon(frames, strand, cdsStart, cdsEnd, exonStart, exonEnd)
    
    
    first_utr_end = cdsStart
    last_utr_start = cdsEnd
    
    if strand == "+" and codon_complete(last_codon):
        cdsEnd = move_pos(txStart, txEnd, exonCount, exonStart, exonEnd, last_utr_start, -3)
    if strand == "-" and codon_complete(first_codon):
        cdsStart = move_pos(txStart, txEnd,exonCount, exonStart, exonEnd, first_utr_end, 3)
        

    if geneName not in GENES.keys():
        build_gene_line(source, geneName, chrom, strand, "gene", txStart, txEnd, file)
        GENES[geneName] = True

    build_gtf_line(source, name, chrom, strand, geneName, "transcript", txStart, txEnd, -1, -1, exonCount, file)
    
    for x in range(exonCount):
        build_gtf_line(source, name, chrom, strand, geneName, "exon", exonStart[x], exonEnd[x], x, -1, exonCount, file)
        if cdsStart < cdsEnd:
            write_features(x, source, name, chrom, strand, geneName, first_utr_end, cdsStart, cdsEnd, last_utr_start, frames[x], exonStart, exonEnd, exonCount, file)
            
    if strand != "-":
        if codon_complete(first_codon):
            write_codon(source, name, chrom, strand, geneName, "start_codon", first_codon, exonCount, file)
        if codon_complete(last_codon):
            write_codon(source, name, chrom, strand, geneName, "stop_codon",last_codon, exonCount, file)
    else:
        if codon_complete(last_codon):
            write_codon(source, name, chrom, strand, geneName, "start_codon",last_codon, exonCount, file)
        if codon_complete(first_codon):
            write_codon(source, name, chrom, strand, geneName, "stop_codon", first_codon, exonCount, file)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-bed", 
        help="Path to BED file", 
        required=True,
        type=str
        )
    parser.add_argument(
        "-gtf", 
        help="Path to output GTF file", 
        required=True,
        type=str
        )
    parser.add_argument(
        "-iso", 
        help="Path to isoforms file", 
        required=True,
        type=str
        )

    if len(sys.argv) < 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    return args



def main():
    args = parse_args()
    bed = args.bed
    gtf = open(args.gtf, "w")
    isoforms = hashmap(args.iso)

    start = time.time()
    with open(bed, "r") as file:
        for line in file:
            run(line, gtf, isoforms)
    end = time.time()
    print(f"Time elapsed: {end - start}")



if __name__ == "__main__":
    main()