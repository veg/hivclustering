import csv, argparse, sys, os.path
import hypy as hy

#-------------------------------------------------------------------------------


arguments = argparse.ArgumentParser(description='Read filenames.')

arguments.add_argument('-i', '--input',   help = 'Input CSV file with cluster definitions. Must be a CSV file with two columns: ID1,cluster ID (integer).', required = True)
arguments.add_argument('-t', '--threshold', help = 'Only report clusters with at least this many sequences')
arguments.add_argument('-a', '--alignment',   help = 'An alignment file with corresponding sequences (must be readable by HyPhy)')
arguments.add_argument('-x', '--extract',   help = 'Extract information about the cluster with a given numeric ID')
arguments.add_argument('-o', '--output',   help = 'Output sequences grouped by cluster to this directory, naming them seqid.nex (default is to output to the current directory)')

settings = arguments.parse_args()

try:
	settings.input = open (settings.input, 'r')
except IOError:
	print ("Failed to open '%s' for reading" % (settings.input))
	raise

	
if settings.alignment != None:
	try:
		test = open (settings.alignment, 'r')
		test.close() 
	except IOError:
		print ("Failed to open '%s' for reading" % (settings.alignment))
		raise
		
if settings.output != None:
	if os.path.isdir (settings.output) == False:
		raise Exception ("'%s' is not an existing directory " % (settings.output))
else:
	settings.output = os.getcwd()


if settings.threshold == None:
	settings.threshold = 1
else:
	settings.threshold = int (settings.threshold)
	
if settings.extract != None:
	settings.extract = int (settings.extract)

clusterReader = csv.reader(settings.input)
header = next(clusterReader)

if len (header) != 2:
	raise Exception ('Expected a .csv file with 2 columns as input')
	
clusters = {}

for line in clusterReader:
	clusterID = int(line[1])
	if settings.extract == None or clusterID == settings.extract:
		if clusterID not in clusters:
			clusters[clusterID] = [line[0],]
		else:
			clusters[clusterID].append (line[0])
		
write_these = {}	
		
for cluster in clusters:
	if len (clusters[cluster]) >= settings.threshold:
		print ("\nCluster %d with %d sequences" % (cluster, len(clusters[cluster])))
		write_these [cluster] = ";".join(clusters[cluster])
		
if settings.alignment != None:
	print ("Processing the alignment through HyPhy")
	hy_instance = hy.HyphyInterface ();
	
	script_path = os.path.realpath(__file__)
	hbl_path =  os.path.join(os.path.dirname(script_path), "..", "lib", "hivclustering", "data", "HBL", "ExtractACluster.bf")
	
	hy_instance.queuevar ('alignment',  os.path.realpath(settings.alignment))
	hy_instance.queuevar ('cluster_definitions',  write_these)
	hy_instance.queuevar ('output_path', os.path.realpath(settings.output))
	hy_instance.runqueue (batchfile = hbl_path)
	print (hy_instance.stdout)

	
