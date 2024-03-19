#!/usr/bin/env python3
import argparse
import cyvcf2

def arg_parser():
	parser = argparse.ArgumentParser(description='Parse APOE genotypes from VCF file')
	parser.add_argument("-v","--vcf",help="Phased VCF with APOE snps rs429358 and rs7412",required=True)
	parser.add_argument("-g","--genome",choices=["GRCh37","GRCh38"],help="GRCh37 or GRCh38 reference genome",required=True)
	parser.add_argument("-o","--out",help="APOE output file prefix",required=True)
	parser.add_argument("-t","--threads",type=int,help="Number of threads to use",default=1)
	return parser.parse_args()

def get_GT(vcf_reader: cyvcf2.VCF, chr: str, pos: str, samples: list) -> dict:
	gt_list = []
	for record in vcf_reader(chr + ':' + pos + '-' + pos):
		gt_list = record.gt_bases
	if gt_list.size == 0:
		raise ValueError(f'Required SNP {pos} is not found in vcf')
	gt_dict = {}
	for i, sample in enumerate(samples):
		gt_dict[sample] = gt_list[i]
	return gt_dict

def is_phased(gt: dict):
	for sample, genotype in gt.items():
		if '|' not in genotype:
			raise AssertionError(f'Genotype {genotype} for sample {sample} is not phased')

def get_haplotypes(rs429358: dict, rs7412: dict) -> dict:
	haps = {}
	for sample, genotype_429358 in rs429358.items():
		genotype_7412 = rs7412[sample]
		haps[sample] = [
			genotype_429358.split('|')[0] + genotype_7412.split('|')[0],
			genotype_429358.split('|')[1] + genotype_7412.split('|')[1]
		]
	return haps

def apoe_genotyper(haps: dict) -> dict:
	apoe_dict = {'CT': '1', 'TT': '2', 'TC': '3', 'CC': '4'}
	apoe_genotypes = {}
	for sample, haplotypes in haps.items():
		apoe_genotypes[sample] = [
			apoe_dict[haplotypes[0]],
			apoe_dict[haplotypes[1]]
		]
	return apoe_genotypes

def add_e_coding(apoe_genotypes: dict) -> dict:
  e_coding = {}
  for sample, genotypes in apoe_genotypes.items():
    sorted_genotypes = sorted(genotypes)
    e_coding[sample] = [genotypes[0], genotypes[1], 'e' + ''.join(sorted_genotypes)]
  return e_coding

def write_output(apoe_genotypes: dict, out: str):
	with open(out + '.tsv', 'w') as f:
		f.write("ID\tAPOE.0\tAPOE.1\tAPOE.GT\n")
		for sample, genotypes in apoe_genotypes.items():
			f.write(f'{sample}\t{genotypes[0]}\t{genotypes[1]}\t{genotypes[2]}\n')

def main():
	args = arg_parser()
	genome = args.genome
	if genome == 'GRCh37':
		chr, rs429358_pos, rs7412_pos = '19', '45411941', '45412079'
	elif genome == 'GRCh38':
		chr, rs429358_pos, rs7412_pos = 'chr19', '44908684', '44908822'
	vcf_reader = cyvcf2.VCF(fname=args.vcf, lazy=True, threads=args.threads)
	samples = vcf_reader.samples
	rs429358 = get_GT(vcf_reader, chr, rs429358_pos, samples)
	rs7412 = get_GT(vcf_reader, chr, rs7412_pos, samples)
	is_phased(rs429358)
	is_phased(rs7412)
	haps = get_haplotypes(rs429358, rs7412)
	apoe_genotypes = apoe_genotyper(haps)
	e_coded_apoe_genotypes = add_e_coding(apoe_genotypes)
	write_output(e_coded_apoe_genotypes, args.out)

if __name__ == '__main__':
	main()
