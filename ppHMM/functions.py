
def analyse(object_array,filename, profilesDict, logFile, outputDir):
    """Take the csv file and add each line as an ISOLATE object
    into object_array
    """
    from classex1 import ISOLATE
    from collections import defaultdict
    import pygal
    from math import log
    import os
    import numpy as np

    leftDict={}
    for i in profilesDict:
        if len(profilesDict[i])>2:
            leftDict[profilesDict[i][0]]=[i, 'profile1']
            leftDict[profilesDict[i][1]]=[i, 'profile2']

    hmmerOutput = open(filename, 'r')

    E_counter=0
    good_counter=0
    total_counter=0
    e_values = defaultdict(list)  # list of e_values
    counter=0 # counter to miss first line
     #everything in table matched or un matched with e_value over 0.00005
    for item in hmmerOutput :
        counter+=1
        parts=item.split(",")
        total_counter+=1
        profileName=parts[0].replace('\n','')
        if counter >1 and profileName in leftDict.keys():
            e_values[profileName].append(-log(float(parts[11])))
        if counter>1 and float(parts[11])<float(0.00005) and profileName in leftDict.keys():
            good_counter  +=  1
            y = ISOLATE()
            y.found = []
            id_parts = parts[3].split('_')
            y.identifier = '_'.join(id_parts[:2])
            y.HNHstart = int(parts[17])
            y.HNHstop = int(parts[18])
            y.HNHcontig = parts[3]
            y.coltype = leftDict[profileName][0]
            y.found.append(parts[0].replace("\n",""))
            y.contiglength = int(parts[5])
            y.MATCH = leftDict[profileName][1]
            y.found.append(parts[0].replace("\n",""))
            object_array.append(y)
        else:
            E_counter+=1




    # 'Total number of hits: {}<br>Hits that did not pass the E counter threshold: {}<br> \
    #     Passed hits:  {}<br> Unaccounted for hits: {}'.format(total_counter,E_counter,
    #         good_counter,int(total_counter-(E_counter+good_counter))), logFile

    hmmerOutput.close()
    histDict=defaultdict(list)

    #Plotting the e-values as a histogram. Needs a bit of reformatting to deal
    #with pygals plotting format
    for profile in e_values:
        if len(e_values[profile]) < 100:
            binsN = len(e_values[profile])
        else:
            binsN =100
        x = np.histogram(e_values[profile], bins=binsN)
        for i in range(len(x[0])):
            histDict[profile].append((x[0][i], x[1][i], x[1][i+1]))
    hist = pygal.Histogram(stroke=False, x_title= 'E value (-log)')
    hist.title = 'E-value distribution'
    for profile in histDict.keys():
        hist.add(profile, histDict[profile])

    EvaluePath = os.path.join(outputDir, 'E_valueHist.svg')
    hist.render_to_file(EvaluePath)

    logFile.EvalueString(total_counter = total_counter,
                            E_counter = E_counter,
                            good_counter = good_counter,
                            Evaluepygal = EvaluePath)

    return object_array

###############################################################################

###############################################################################

def analyse_dist(object_array,profilesDict, logFile, args):

    if args.single:
        print 'Skipping Distance check as profile is Single'
        matched_array2  =  profilesClass()
        for i in object_array:
            matched_array2.append(i)


    else:
        #Loop through and apply district metric
        from collections import defaultdict
        import sys
        matched_array2  =  profilesClass()
        object_dict=defaultdict(list)
        for i in object_array:
            object_dict[i.HNHcontig[:-1]].append(i)
        counter=0
        for K in object_dict.keys():

            counter+=1

            for i in object_dict[K]:
                for j in object_dict[K]:
                    if i.coltype==j.coltype:
                        if i.MATCH=='profile2' and j.MATCH=='profile1':

                            if i.HNHstart-j.HNHstop >=0 and i.HNHstart-j.HNHstop<=int(profilesDict[i.coltype][2]):
                                j.IMMstart  = i.HNHstart
                                j.IMMstop   = i.HNHstop
                                j.IMMcontig = i.HNHcontig

                                j.MATCH="TRUE"

                                matched_array2.append(j)

        profiles_count = matched_array2.count_profiles()
    return  matched_array2






def plot_dist(object_array, profilesDict, outputDir, logFile, args):
    from collections import defaultdict
    import sys
    import os
    import pygal
    from pygal.style import LightStyle
    matched_array2=[]
    object_dict=defaultdict(list)
    for i in object_array:
        object_dict[i.HNHcontig[:-1]].append(i)
    counter=0
    countingDict=defaultdict(dict)
    plottingDict=defaultdict(list)
    for K in object_dict.keys():

        counter+=1

        for i in object_dict[K]:
            for j in object_dict[K]:
                if i.coltype==j.coltype:
                    if i.MATCH=='profile2' and j.MATCH=='profile1' and i.HNHstart-j.HNHstop > -int(profilesDict[i.coltype][2]) \
                    and i.HNHstart-j.HNHstop < 3*int(profilesDict[i.coltype][2]) :
                        try:
                            countingDict[i.coltype][i.HNHstart-j.HNHstop]+=1
                        except:
                            countingDict[i.coltype][i.HNHstart-j.HNHstop]=0

    for i in countingDict:
        for j in countingDict[i]:
            plottingDict[i].append((countingDict[i][j], j,j+1))
    hist = pygal.Histogram()
    hist.title = 'Hits found across inter-genic region'
    for i in plottingDict:
        hist.add(i, plottingDict[i])

    interProfilePath = os.path.join(outputDir, 'InterProfileHist.svg')
    hist.render_to_file(interProfilePath)

    logFile.interProfile(args, interProfilePath)





def choose_profiles(matched_array):
    #given an arry of hits this function will set up a GUI to select the allowed profiles and return them as strings in a vector
    old_profiles=[]
    for i in matched_array:
        for j in i.profiles:
            if j not in old_profiles:
                old_profiles.append(j)
    global NewProfile_list
    NewProfile_list=[]

    label_dict={}
    import Tkinter as tk
    master=tk.Tk()

    def ammend_list():
        for i in label_dict.keys():
            if label_dict[i].get()==1:
               NewProfile_list.append(i)


    tk.Label(master,
    text="Profile selection! Select profiles you want removing!").grid(row=0, sticky=tk.W)
    for i in range(len(old_profiles)):
        label_dict[old_profiles[i]]=tk.IntVar()
        j = (i/8)
        k=i-(8*j)
        tk.Checkbutton(master, text=old_profiles[i], variable=label_dict[old_profiles[i]]).grid(row=k+1, column=j, sticky=tk.W)

    #def ammend_list():
    #    for i in old_profiles:
    tk.Button(master, text="Select", command=ammend_list).grid(row=i+2, stick=tk.W)
    tk.Button(master, text='Quit', command=master.quit).grid(row=i+3, sticky=tk.W, pady=4)

    master.mainloop()
    profile_list=NewProfile_list
    return profile_list

###############################################################################

###############################################################################

class profilesClass(list):

    def __init__(self):
        self.test='TEST'


    def make_addprofiles(self):
        file_handle=open('Seqs_Pfam',"w")

        for i in self:
            file_handle.write('>{}\n{}\n'.format(i.seqid, i.HNHsequence))
        file_handle.close()

    def search_profiles(self, args):
        import os
        import subprocess
        if args.remove:
            PfamDir = args.remove
        elif args.cluster:
            PfamDir = args.cluster

        print 'Searching sequences against {}'.format(PfamDir)

        command ='hmmscan --tblout Pfam_annotation {} Seqs_Pfam'.format(PfamDir)
        process = subprocess.Popen(command.split(' '),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        process.communicate()


    def addprofiles(self):
        from collections import defaultdict
        file_handle=open('Pfam_annotation',"r")
        count=0
        annoDict=defaultdict(list)
        for line in file_handle:

            if line.startswith('#')==False:
                    parts= [x for x in line.split(' ') if x!='']

                    if float(parts[4])<float(0.005):

                        annoDict[int(parts[2])].append(parts[0])
            for i in self:
                try:
                    i.profiles = annoDict[i.seqid]
                except:
                    pass
        file_handle.close()

    def seqID(self):
        counter=0
        for i in self:
            counter+=1
            i.seqid=counter

    def clusterSeqs(self):
        import subprocess
        #run cd_hit command on the sequences
        file_handle=open('Cluster_seqs.fa','w')
        for i in self:
            string=">{}\n{}\n".format(i.seqid, i.HNHsequence)
            file_handle.write(string)
        file_handle.close()

        search = "cd-hit -i Cluster_seqs.fa -o Clustered -c 0.98 -d 0"
        command = search.split(' ')
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        process.communicate()
        #Read the cd-hit file and add sequences to the correct cluster
        Cluster_dict = {}
        Lead_cluster = {}
        file_handle = open("Clustered.clstr",'r')
        for line in file_handle:
            if line.startswith('>'):
                Cluster = ''.join(line[9:])
                Cluster = Cluster.strip('\n')
            else:
                line = line.replace('.','>')
                parts = line.split('>')
                Cluster_dict[int(parts[1])]=Cluster
                if '*' in line:
                    Lead_cluster[int(parts[1])]='Lead'
                else:
                    Lead_cluster[int(parts[1])]='follow'


        #Make one sequence a 'Cluster leader'

        for i in self:
            i.Cluster = Cluster_dict[i.seqid]
            i.Lead    = Lead_cluster[i.seqid]

    def plot_profiles(self, filename, logFile):
        from collections import defaultdict
        import pygal
        line_chart = pygal.Bar()
        line_chart.title = 'Pfam profiles associated with predicted sequences'
        profile_dict = defaultdict(dict)
        profile_list = []
        for i in self:
            for j in i.profiles:
                try:
                    profile_dict[i.coltype][j]+=1
                except:
                    profile_dict[i.coltype][j]=1
                    profile_list.append(j)
                    pass
        profile_list = list(set(profile_list))
        line_chart.x_labels = profile_list
        for i in profile_dict:
            plotList=[]
            for k in profile_list:
                try:
                    plotList.append(profile_dict[i][k])
                except:
                    plotList.append(0)
            line_chart.add(i, plotList)

        line_chart.render_to_file(filename)
        logFile.profilePygal = filename

    def finalOutPut(self,outputDir, args):
        #Write the results in a csv file with sequences
        import os
        Passed_handle = open(os.path.join(outputDir,'Passed_output.csv'),'w')
        if args.remove:
            header = 'Identifier,Pass/Fail,Species,Isolate,Bacteriocin,Immunity,Profiles \n'
            Passed_handle.write(header)
            for i in self:
                Passed_handle.write(','.join([i.identifier, str(i.profilePass), i.Species, i.isolate
                ,str(i.HNHsequence), str(i.IMMsequence),str(i.profiles).replace('[','').replace(']','')]))
                Passed_handle.write('\n')
        elif args.cluster:
            header = 'Identifier,Pass/Fail,Species,Isolate,Cluster,Bacteriocin,Immunity,Profiles \n'
            Passed_handle.write(header)
            for i in self:
                Passed_handle.write(','.join([i.identifier,str(i.profilePass), i.Species, i.isolate
                ,i.Cluster, str(i.HNHsequence), str(i.IMMsequence),str(i.profiles).replace('[','').replace(']','')]))
                Passed_handle.write('\n')
        else:
            header = 'Identifier,Species,Isolate,Bacteriocin,Immunity \n'
            Passed_handle.write(header)
            for i in self:
                Passed_handle.write(','.join([i.identifier, i.Species, i.isolate
                ,str(i.HNHsequence), str(i.IMMsequence)]))
                Passed_handle.write('\n')
        Passed_handle.close()

        #Now make a couple of fasta files with Sequences
        PassedFastaN = open(os.path.join(outputDir,'PassedFastaN.fa'), 'w')
        PassedFastaC = open(os.path.join(outputDir,'PassedFastaC.fa'), 'w')
        FailedFastaN = open(os.path.join(outputDir,'FailedFastaN.fa'), 'w')
        FailedFastaC = open(os.path.join(outputDir,'FailedFastaC.fa'), 'w')
        for i in self:
            fastaStringN = str('>'+i.identifier+'_'+str(i.seqid)+'\n'+i.HNHsequence+'\n')
            fastaStringC = str('>'+i.identifier+'_'+str(i.seqid)+'\n'+i.IMMsequence+'\n')
            if i.profilePass:
                PassedFastaN.write(fastaStringN)
                PassedFastaC.write(fastaStringC)
            else:
                FailedFastaN.write(fastaStringN)
                FailedFastaC.write(fastaStringC)

    def seqLengthPlot(self, outputDir,logFile):
        from collections import defaultdict
        import numpy as np
        import pygal
        import os
        #Plotting the e-values as a histogram. Needs a bit of reformatting to deal
        #with pygals plotting format
        e_values = defaultdict(list)
        histDict = defaultdict(list)
        for i in self:
            e_values[i.coltype].append(len(i.HNHsequence))
        for profile in e_values:
            if len(e_values[profile]) < 100:
                binsN = len(e_values[profile])
            else:
                binsN =100
            x = np.histogram(e_values[profile], bins=binsN)
            for i in range(len(x[0])):
                histDict[profile].append((x[0][i], x[1][i], x[1][i+1]))
        hist = pygal.Histogram(stroke=False, x_title= 'Sequence length (aa)')
        hist.title = 'Sequence length distribution'
        for profile in histDict.keys():
            hist.add(profile, histDict[profile])

        EvaluePath = os.path.join(outputDir, 'seqLengthHist.svg')
        hist.render_to_file(EvaluePath)
        logFile.seqDist(EvaluePath)






    def seqFinder_prodigal(self, args):
        #Will find the corordinates in translation and find the correct open
        #reading frame
        import os
        from Bio import SeqIO
        import subprocess

        genome_dir    =  args.Genome_dir
        orf_dir  =  args.ORF_dir
        counter  =  0
        print 'Finding sequences...'
        print len(self)

        for search in self:
            counter+=1

            isolate_search=search.identifier

            if args.METAfile:
                pass
            else:
                search.contig_file = str(search.identifier)

            #start and end points of motifs
            endposHNH  =  search.HNHstop
            startposHNH  =  search.HNHstart


            #get the contigs
            HNHcontig  =  search.HNHcontig
            collect  =  ""
            #FRAME
            frameHNH=int(HNHcontig[-1:])

            #open up the sequence file
            parts=HNHcontig.split("_")
            genome_file  =  '{}{}'.format(genome_dir, search.contig_file)
            ORF_file  =  '{}{}'.format(orf_dir,HNHcontig.replace('\n',''))

            if os.path.exists(genome_file)==False:
                print "File doesnt exit!!"
                print genome_file
                search.finder = "NO_FILE"
                continue
            try:
                file_handle=open(genome_file,"r")
            except:
                continue



            #open the genome and find the contif file
            output_file  =  'seqtest.fa'
            records  =  (r for r in SeqIO.parse(genome_file, "fasta") if HNHcontig[:-2] in r.id.upper())
            count  =  SeqIO.write(records, output_file, "fasta")

             # runs transeq on the command line
            command  =  "transeq -sequence {} -outseq _6frame.fa -clean -table 11 -frame 6".format(genome_file)
            commandTranseq = command.split(' ')

            #run prodigal ORF finding software on contig
            command  = ("prodigal -i seqtest.fa -a {} -d tester2.fa -f gbk -o output_file_waste.txt"
                " -s output_prod_genes.txt -g 11 -p meta").format(ORF_file)
            command_getorf = command.split(' ')


            #for now run all files again
            if os.path.isfile(ORF_file):
                replace_spaces  =  'sed -i -e "s/ //g" ' + (ORF_file)
                replace_spaces = replace_spaces.split(' ')
                p = subprocess.Popen(replace_spaces, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                p.communicate()
            else:
                p = subprocess.Popen(command_getorf, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                p.communicate()

                replace_spaces  =  'sed -i -e "s/ //g" {}'.format(ORF_file)
                replace_spaces = replace_spaces.split(' ')
                p = subprocess.Popen(replace_spaces, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                p.communicate()

                replace_spaces  =  'sed -i -e "s/ //g" tester2.fa'
                replace_spaces = replace_spaces.split(' ')
                p = subprocess.Popen(replace_spaces, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                p.communicate()

            p = subprocess.Popen(commandTranseq, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            p.communicate()
            #Find the correct translation from the transeq and add it to string
            input_file  =  '_6frame.fa'
            records = (r for r in SeqIO.parse(input_file, "fasta") if HNHcontig in r.id.upper())
            for i in records:
                stringHNH  =  i.seq

            #find the colicin and the genes either side if in single mode or the immunity if in multi mode
            search.IMMsequence_1=""
            search.error='No_error'
            try:
                records = list(SeqIO.parse(ORF_file, "fasta"))
            except:
                print 'Error in making ORF file for {}'.format(search.identifier)
                continue
            hmmerString = str(stringHNH[startposHNH+5:endposHNH-5])
            for i in range(len(records)):

                if hmmerString in records[i].seq:
                    search.HNHsequence  =  records[i].seq
                    taken_id = records[i].id

                    try:
                        parts  =  records[i].id.split('#')
                        strand  =  parts[3]
                        search.strand  =  parts[3]
                        if strand=='1':
                            C_take   =  1
                            N_take   =  -1
                        elif strand  ==  '-1':
                            N_take   =  1
                            C_take   =  -1

                        search.HNHsequence_prokkaend  =  parts[2]
                        search.HNHsequence_prokkastart  =  parts[1]


                        parts_C  =  records[i+C_take].id.split('#')
                        if parts_C[3]  ==  strand:

                            search.IMMsequence  =  records[i+C_take].seq
                            parts  =  records[i+C_take].id.split('#')
                            search.IMMsequence_prokkastart  =  parts[1]

                        else:

                            search.IMMsequence  =  records[i+(2*C_take)].seq
                            parts  =  records[i+(2*C_take)].id.split('#')
                            search.IMMsequence_prokkastart  =  parts[1]

                        search.IMMsequence_1  =  records[i+N_take].seq
                        parts  =  records[i+N_take].id.split('#')
                        search.IMMsequence_1_prokkastart  =  parts[1]
                        search.error  =  'No_error'
                        print 'Found {}'.format(search.identifier)
                    except:
                        print 'Error in {}'.format(search.identifier)
                        search.error  =  "CONTIG_ENDED"
                        pass



    def count_profiles(self):
        from collections import defaultdict
        profile_count  =  defaultdict(int)
        for i in self:
            profile_count[i.coltype]+=1

        return profile_count
