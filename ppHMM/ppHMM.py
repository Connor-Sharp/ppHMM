def main():

	from argparse import ArgumentParser
	import cPickle as pickle
	from html_writer import html_writer
	import os
	from time import gmtime, strftime
	from collections import defaultdict

	#get current directory and make an output folder (including html)
	cwd = os.getcwd()
	date = strftime("%Y_%m_%d_%H_%M", gmtime())
	outputDir = os.path.join(cwd, date)
	outputDir+='_'
	Counter=0
	outputDir+=str(Counter)
	while os.path.isdir(outputDir):
		outputDir = outputDir[:-1]
		outputDir+=str(Counter)
		Counter+=1
	os.mkdir(outputDir)

	logFileHandle = open(os.path.join(outputDir, 'Summary.html'), 'w')
	logFile = html_writer(logFileHandle)

	parser = ArgumentParser(description='This script will identify the PFAM profiles in the hmmer output file')

	logFile.writeTitle()
	logFile.makeTabs()
	#Define optional and required arguments
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	# Add more options if you like
	required.add_argument("-fin", "--file_in", dest="myFilenameIn",
	                    help="File from HMMer output", metavar="FILE", required=True)

	required.add_argument("-fout", "--file_out", dest="myFilenameOut",
	                    help="Parsed output of hmmerfile", metavar="FILE", required=True)

	required.add_argument("-g", "--genomes", dest="Genome_dir",
	                    help="directory of genomes used in hmmer analysis (fasta)", metavar="FILE", required=True)

	required.add_argument("-o", "--orf", dest="ORF_dir",
	                    help="Directory to place ORF predictions from fasta genomes", metavar="FILE", required=True)

	optional.add_argument("-m", "--meta", dest="METAfile",
	                    help="File with id, ISOLATE, Species and contigFile", metavar="FILE")

	required.add_argument("-p", "--profile", dest="profile_file", help="Text file with names of profiles to be identified formatted as: \
		NAME,1stprofile,2ndprofile,interprofiledistance", required=True)

	optional.add_argument("-r","--remove", dest="remove", metavar="FILE",
		help='Remove sequences based on unwanted Pfam domain', required=False)

	optional.add_argument("-rc","--Cluster", dest="cluster", metavar="FILE",
		help='Sequences will be clustered (0.95 seqid) before pfam removal', required=False)
	optional.add_argument("-s","--Csingle", dest="Single", metavar="FILE",
		help='Sequences will be clustered (0.95 seqid) before pfam removal', action='store_true',required=False)

	args = parser.parse_args()

	logFile.writeArgs(args)

	def read_profiles(profile_file):
		#Read through the profile file which has all the NAMES and profiles in a comma seperated file
		profileDict={}
		file_handle=open(profile_file)
		for line in file_handle:
			parts = line.strip('\n').split(',')
			profileDict[parts[0]] = parts[1:]
		return profileDict



	def read_hmmerinput(hmmerInput,hmmerOutput, log):
		#Read the space seperated output of hmmer and change it into a much more usable csv table.

		file_handle = open(hmmerInput,'r')

		file_output = open(hmmerOutput,'w')
		file_output.write(str("Target_name,accession,tlen,query_name,species,qlen,E-value,score,\
	bias,#,of,cEvalue,iE-value,score,bias,from,to,from,to,from,to,acc,Description\n"))
		for line in file_handle:
			if line.startswith('#')==False:
						 	parts=line.split(' ')
						 	parts=[x for x in parts if x !='']
							file_output.write(','.join(parts[:22]))
							file_output.write(',')
							file_output.write(' '.join(parts[22:]))
		file_handle.close()
		file_output.close()


	#read profiles from file
	profileDict={}
	profileDict = read_profiles(args.profile_file)


	logFile.profileTableWrite(profileDict)


	read_hmmerinput(args.myFilenameIn, args.myFilenameOut, logFile)


	from functions import analyse
	object_array=[]
	object_array = analyse(object_array, args.myFilenameOut, profileDict, logFile, outputDir)

	from functions import plot_dist
	plot_dist(object_array,profileDict,outputDir,logFile)

	#from functions import addBIGSspecies
	#addBIGSspecies(object_array, args, logger)


	from functions import profilesClass
	matched_array=profilesClass()
	temp_array=profilesClass()


	from functions import analyse_dist
	matched_array = analyse_dist(object_array,profileDict, logFile)
	test_arrayClean = profilesClass()
	#id,isolate,species ,Genome file
	if args.METAfile:
		metaDict=defaultdict(list)
		file_handle = open(args.METAfile,'r')
		for line in file_handle:
			parts = line.split(',')
			metaDict[parts[0]] = parts[1:]
		for i in matched_array:
			try:
				infoList = metaDict[i.identifier]
				i.Species = infoList[1]
				i.isolate = ''.join([infoList[1], infoList[0]])
				i.contig_file = infoList[2].replace('\n','')
			except:
				print 'Missing data for {}'.format(i.identifier)






	matched_array.seqFinder_prodigal(args)

	for i in matched_array:
		if len(i.HNHsequence) > 10:
			test_arrayClean.append(i)

	matched_array = test_arrayClean
	matched_array.seqID()







	#pickle_file_binary=str("results_binary_profiles.pickle")
	#pickle_file=str("results_profiles.pickle")

	#with open(pickle_file_binary, "wb") as output_file:
	#    pickle.dump(match_it, output_file)

	#with open(pickle_file, "w") as output_file:
	#    pickle.dump(match_it, output_file)

	if args.remove or args.cluster:
		from functions import choose_profiles
		matched_array.make_addprofiles()
		matched_array.search_profiles(args)
		matched_array.addprofiles()
		matched_array.plot_profiles(os.path.join(outputDir, 'All_profiles.svg'), logFile)
		profile_list = choose_profiles(matched_array)
		passed, failed = [0,0]
		if args.cluster:
			matched_array.clusterSeqs()

		if args.remove:
			#remove sequences based on profiles
			for i in matched_array:
				i.profilePass = True
				for j in i.profiles:
					if j in profile_list:
						i.profilePass = False



		if args.cluster:
			ClusterPass={}
			#remove clusters by profile list
			for i in matched_array:
				if i.Lead =='Lead':
					ClusterPass[i.Cluster] = True
					for j in i.profiles:
						if j in profile_list:
							ClusterPass[i.Cluster] = False

			for i in matched_array:
				i.profilePass = ClusterPass[i.Cluster]


			for i in matched_array:
				if i.profilePass:
					passed +=1
				else:
					failed +=1
				print i.identifier, i.seqid, i.profiles,len(i.HNHsequence),i.profilePass
		logFile.profileAnalysis( args, passed, failed)


	#Plot sequence lengths Passed or failed???



	#Output a csv of all sequences with some information
	matched_array.seqLengthPlot(outputDir, logFile)
	matched_array.finalOutPut(outputDir,args)

	#Output a pickle file of the profileClass() object for future use.
	pickle_file_Passed=os.path.join(outputDir,"resultsPassed.pickle")
	pickle_file_Failed=os.path.join(outputDir,"resultsFailed.pickle")

	matchPassed =[]
	matchPassed = [x for x in matched_array if x.profilePass]
	matchFailed = []
	matchFailed = [x for x in matched_array if x.profilePass == False]
	with open(pickle_file_Passed, "w") as output_file:
	   pickle.dump(matchPassed, output_file)

	with open(pickle_file_Failed, "w") as output_file:
	   pickle.dump(matchFailed, output_file)

	logFile.finish()

if __name__ == '__main__':
    main()
