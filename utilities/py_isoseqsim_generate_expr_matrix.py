#!/usr/bin/env python
import sys,time,argparse
import numpy as np
from numpy.random import multinomial

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_expr_matrix(args.input,args.expression,args.nb_n,args.nb_p,args.output,args.sum_counts)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
def resample(total_num_reads,read_counts):
	return multinomial(total_num_reads, read_counts/sum(read_counts))
def generate_expr_matrix(input_gpd_fl,expression_fl,nb_r,nb_p,output_expr_mtx,total_num_reads):
	# parse gpd file
	gpd_dict = {}
	for line in input_gpd_fl:
		isoform_name = line.split('\t')[1].strip()
		gpd_dict[isoform_name] = line.strip()
	read_count_list = []
	isoform_list = []
	for line in expression_fl:
		[isoform_name,count] = line.split('\t')[0].strip(),line.split('\t')[1].strip()
		read_count_list.append(float(count))
		isoform_list.append(gpd_dict[isoform_name])
	resampled_read_count_list = resample(total_num_reads,np.array(read_count_list)).tolist()
	for read_count,line in zip(resampled_read_count_list,isoform_list):
		print >>output_expr_mtx, line + "\t" + str(read_count)
	# # generate random read count based on negative binomial distribution
	# nb_list = np.random.negative_binomial(nb_r,nb_p,len(gpd_list)).tolist()
	# i = 0
	# for gpd in gpd_list:
	# 	print >>output_expr_mtx, gpd + "\t" + str(nb_list[i])
	# 	i += 1
	input_gpd_fl.close()
	expression_fl.close()
	output_expr_mtx.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Randomly generate read count for each isoform based on negative binomial (NB) distribution. Read count is shown in last column of output file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-expr','--expression',type=argparse.FileType('r'),required=True,help="Input: expression file")
	parser.add_argument('-n','--nb_n',type=int,default=10,help="Parameter of the NB distribution, n")
	parser.add_argument('-p','--nb_p',type=float,default=0.5,help="Parameter of the NB distribution, p")
	parser.add_argument('--sum_counts',type=int,help='Total number of read count')
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd + read count file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
