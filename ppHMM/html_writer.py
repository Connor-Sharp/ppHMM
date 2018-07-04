

class html_writer:

	header = """
	<!DOCTYPE html>
	<html>
	<head>
	  <meta charset="utf-8">
	  <meta name="viewport" content="width=device-width, initial-scale=1">
	  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
	  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
	  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
	<title>Output</title>

	<style>
	table {
	    font-family: arial, sans-serif;
	    border-collapse: collapse;
	    width: 40%;
	}

	td, th {
	    border: 1px solid #dddddd;
	    text-align: center;
	    padding: 4px;
	}

	tr:nth-child(even) {
	    background-color: #dddddd;
	}

	p {
	    margin: 35px;
	}
	/* Style the tab */
.tab {
    overflow: hidden;
    border: 1px solid #ccc;
    background-color: #f1f1f1;
}

/* Style the buttons inside the tab */
.tab button {
    background-color: inherit;
    float: left;
    border: none;
    outline: none;
    cursor: pointer;
    padding: 14px 16px;
    transition: 0.3s;
    font-size: 17px;
	}

	/* Change background color of buttons on hover */
	.tab button:hover {
	    background-color: #ddd;
	}

	/* Create an active/current tablink class */
	.tab button.active {
	    background-color: #ccc;
	}

	/* Style the tab content */
	.tabcontent {
	    display: none;
	    padding: 6px 12px;
	    border: 1px solid #ccc;
	    border-top: none;
	}
	</style>
	</head>
	<body>
	"""
	Tail = """
		<script>
		function openTab(evt, cityName) {
		    var i, tabcontent, tablinks;
		    tabcontent = document.getElementsByClassName("tabcontent");
		    for (i = 0; i < tabcontent.length; i++) {
		        tabcontent[i].style.display = "none";
		    }
		    tablinks = document.getElementsByClassName("tablinks");
		    for (i = 0; i < tablinks.length; i++) {
		        tablinks[i].className = tablinks[i].className.replace(" active", "");
		    }
		    document.getElementById(cityName).style.display = "block";
		    evt.currentTarget.className += " active";
			}
			// Get the element with id="defaultOpen" and click on it
			document.getElementById("defaultOpen").click();
		</script>
		</body>
		</html>
		"""

	def __init__ (self, file_handle=''):
		self.file_handle = file_handle
		self.file_handle.write(self.header)
		self.profilePygal =''



	def writeTitle(self):
		from time import gmtime, strftime
		title = """
		<div class="jumbotron text-center">
		  <h1>Output of ppHMM results</h1>
		  <p>Script run at: {}</p>
		</div>

		""".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
		self.file_handle.write(title)

	def makeTabs(self):
		tabHeader='''
		<div class="tab">
		  <button class="tablinks" onclick="openTab(event, 'Input')" id="defaultOpen">Input and Output</button>
		  <button class="tablinks" onclick="openTab(event, 'HMMerStatistics')">Hmmer Statistics</button>
		  <button class="tablinks" onclick="openTab(event, 'Interprofile')">Inter-profile distances</button>
		  <button class="tablinks" onclick="openTab(event, 'profileAnalysis')">Profile analysis</button>
		  <button class="tablinks" onclick="openTab(event, 'SeqLength')">Sequence length distribution</button>
		</div>
		'''
		self.file_handle.write(tabHeader)
	def writeArgs(self, args):
		StartTab = """
		<div id='Input' class="tabcontent">
			<h3>Input and Output</h3>
		"""

		hmmerFiles = '<p>Opened HMMer output: {} and converted it into .csv: {}<br>'.format(args.myFilenameIn, args.myFilenameOut)
		Directories = 'Directory of genome files: {}<br>Directory of ORFS: {}<br>'.format(args.Genome_dir, args.ORF_dir)
		if args.cluster:
			pfamFiles = 'Pfam hmmer file: {}<br>'.format(args.cluster)
		elif args.remove:
			pfamFiles = 'Pfam hmmer file: {}<br>'.format(args.remove)
		else:
			pfamFiles = 'Pfam hmmer file: No Pfam file<br>'
		if args.METAfile:
			metaFile = 'Metadata File: {}<br>'.format(args.METAfile)
		else:
			metaFile = 'Metadata File: No Metadata file given<br>'
		self.file_handle.write(StartTab)
		self.file_handle.write(hmmerFiles)
		self.file_handle.write(Directories)
		self.file_handle.write(pfamFiles)
		self.file_handle.write(metaFile)
		self.file_handle.write('</p>')

	def EvalueString(self, total_counter, E_counter, good_counter, Evaluepygal):
		StartTab = """
		<div id='HMMerStatistics' class="tabcontent">
			<h3>Hmmer Statistics</h3>
		"""
		Evalues = '<p>Total number of hits: {}<br>Hits that did not pass the E counter threshold: {}<br>Passed hits:  {}<br> Unaccounted for hits: {}</p>'.format(total_counter,E_counter,good_counter,int(total_counter-(E_counter+good_counter)))
		self.file_handle.write(StartTab)
		self.file_handle.write(Evalues)
		self.pygalWrite(Evaluepygal)
		self.file_handle.write('</div>')



	def writehtmlString(str, file_handle):
		self.file_handle.write("<p>")
		self.file_handle.write(str)
		self.file_handle.write("</p>")

	def pygalWrite(self, imageFile):

		self.file_handle.write('<embed type="image/svg+xml"  src=\"{}\" width=\"900\" height=\"900\">'.format(imageFile))


	def profileTableWrite(self,profileDict):
		tableString = """
	<p>Table of profiles used in this run:</p>

	<table align = "center">
	  <tr>
	    <th>Profile pair</th>
	    <th>N-terminal profile</th>
	    <th>C-terminal profile</th>
	    <th>Inter-profile distance (aa) </th>
	  </tr>


	"""
		self.file_handle.write(tableString)
		outerString =  """
		<tr>
		{}
		</tr>
		"""
		innerString = "<td>{}</td>"

		for Key in profileDict:
			rowList= []
			rowList.append(innerString.format(Key))
			for value in profileDict[Key]:
				rowList.append(innerString.format(value.strip('\n')))
			rowList.append("</th>")
			self.file_handle.write(outerString.format('\n'.join(rowList)))
		self.file_handle.write("</table>")
		self.file_handle.write('</div>')

	def interProfile(self, args, pygalPath):
		StartTab = """
		<div id='Interprofile' class="tabcontent">
			<h3>Inter-profile Distances</h3>
		"""
		self.file_handle.write(StartTab)
		if args.single:
			self.file_handle.write('<p>You selected to scan for only a single profile</p>')
		else:
			self.pygalWrite(pygalPath)
		self.file_handle.write('</div>')

	def finish(self):
		self.file_handle.write(self.Tail)

	def profileAnalysis(self,args, passed, failed):
		StartTab = """
		<div id='profileAnalysis' class="tabcontent">
			<h3>Profiles identified</h3>
		"""
		self.file_handle.write(StartTab)
		if args.cluster or args.remove:
			profileString = 'Sequences which passed profile checks: {}<br>Sequences which failed profile checks: {}'.format(passed, failed)
			self.file_handle.write(profileString)
			self.pygalWrite(self.profilePygal)
		else:
			profileString = 'No profile selection chosen (options -r or -rc)'
			self.file_handle.write(profileString)
		self.file_handle.write('</div>')


	def seqDist(self, SeqpygalPath):
		StartTab = """
		<div id='SeqLength' class="tabcontent">
			<h3>Sequence size distribution</h3>
		"""
		self.file_handle.write(StartTab)
		self.pygalWrite(SeqpygalPath)
		self.file_handle.write('</div>')
