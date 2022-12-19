import os

file_loc = r'/scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_PROW/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals.forbams_final.bed' # change this to your input file location
folder = r'/scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_PROW' # change this to the location you want the output file saved

node_list = {}
node_lengths = {}
max_node = 0

with open(file_loc, 'r') as rfile:
	for line in rfile:
		cols = line.split() # splits entire line of input file where there is a space
		node = cols[0] # column 1: node name
		position = int(cols[1]) # column 2: position
		node_data = node.split('_') # splits the node name into pieces
		node_order = int(node_data[1]) # second piece of node name is its number
		node_length = int(node_data[3]) # fourth piece of node name is its length

		if node_order not in node_list:
			node_list[node_order] = [position] # node_list dictionary keeps track of all positions per node
		else:
			node_list[node_order].append(position)
		
		if node_order not in node_lengths: # node_lengths dictionary keeps track of the length of each node
			node_lengths[node_order] = node_length
		
		if node_order > max_node: # max_node variable is used to set up the for loop in the next step
			max_node = node_order
			
position_list = [] # this list is basically building each row of the output file
length = 0 # variable to keep track of length as we go
			
for node in range(1, max_node + 1): # for loop is iterating through node numbers. same node number used to access both dictionaries
	try:
		positions = node_list[node] # list of positions for that node
		for pos in positions:
			position_list.append(pos + length)
		length += node_lengths[node] # node length is added to the total AFTER that node's positions have been calculated
	except KeyError as e:
		continue
	
			
filename = 'PROW_ref_reduced_concat_snpfile.snp' # change this to your desired output filename
filepath = os.path.join(folder, filename)
with open(filepath, 'w') as wfile:
	for item in position_list:
		wfile.write(str(item))
		wfile.write('\n')
			