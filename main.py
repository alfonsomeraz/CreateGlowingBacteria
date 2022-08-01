# Project: Create Glowing Bacteria

def add_complement(original_DNA_strand):
    # takes DNA strand and returns a complementary DNA strand
    complement_strand = []
    bases = list(original_DNA_strand)
    for base in bases:
        if base == 'A':
            complement_strand.append('T')
        elif base == 'T':
            complement_strand.append('A')
        elif base == 'G':
            complement_strand.append('C')
        elif base == 'C':
            complement_strand.append('G')
    return ''.join(complement_strand)


def cut_plasmid(strand, restriction_site):
    # takes an original and complementary pair of DNA strands and returns
    # them cut at the restriction sites
    original = strand
    complement = add_complement(strand)
    complement_restriction_site = add_complement(restriction_site)

    original_cut = original.find(restriction_site)
    complement_cut = complement.find(complement_restriction_site)

    original_begin = original[:original_cut + 1]
    original_end = original[original_cut + 1:]
    complement_begin = complement[:complement_cut + 5]
    complement_end = complement[complement_cut + 5:]
    return original_begin, original_end, complement_begin, complement_end


def cut_GFP(strand, restriction_sites):
    # takes in GFP strand and restriction sites to create complements
    # cuts strands at restriction sites and returns newly cut GFP strands
    first_restriction = restriction_sites[0]
    second_restriction = restriction_sites[1]
    complement = add_complement(strand)
    first_restriction_comp = add_complement(first_restriction)
    second_restriction_comp = add_complement(second_restriction)

    original_GFP = strand[strand.find(first_restriction)+1:]
    original_GFP = original_GFP[:original_GFP.rfind(second_restriction)+1]

    complement_GFP = complement[complement.find(first_restriction_comp)+5:]
    complement_GFP = complement_GFP[:complement_GFP.rfind(second_restriction_comp)+5]

    return original_GFP, complement_GFP


def ligation(dna_strands, gfp_strands):
    original_first_strand = dna_strands[0]
    original_second_strand = dna_strands[1]
    complement_first_strand = dna_strands[2]
    complement_second_strand = dna_strands[3]

    original = original_first_strand + gfp_strands[0] + original_second_strand
    complement = complement_first_strand + gfp_strands[1] + complement_second_strand
    return original, complement


if __name__ == '__main__':
    filename = input()
    file = open(filename, 'r')
    strands = [i.strip('\n') for i in file.readlines()]
    file.close()

    reads = int(len(strands) / 4)
    data = {}
    new_GFP_strands = []

    for i in range(reads):
        data[i] = {}
        data[i]['original_plasmid_strand'] = strands[(4*i)+0]
        data[i]['restriction_site_plasmid'] = strands[(4 * i) + 1]
        data[i]['original_GFP_strand'] = strands[(4 * i) + 2]
        data[i]['restriction_sites'] = strands[(4 * i) + 3]

    for i in range(len(data)):
        # Get strand details from dictionary
        plasmid = data[i]['original_plasmid_strand']
        plasmid_restriction_site = data[i]['restriction_site_plasmid']
        gfp_strand = data[i]['original_GFP_strand']
        gfp_restriction_sites = data[i]['restriction_sites'].split()

        # cut plasmid and gfp strands
        plasmid_strands = cut_plasmid(plasmid, plasmid_restriction_site)
        gfp_strands = cut_GFP(gfp_strand, gfp_restriction_sites)

        # ligation
        new_strands = ligation(plasmid_strands, gfp_strands)
        new_GFP_strands.append(new_strands)

    for strands in new_GFP_strands:
        print(strands[0])
        print(strands[1])
