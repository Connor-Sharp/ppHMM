# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:32:41 2014

@author: csharp
"""
class ISOLATE:
   'each matched colicin/Im pair is made into an object'
   isolateCount = 0
   def __init__(self,isolate="",coltype="",contiglength=0, identifier="",
                database="",HNHcontig="",IMMcontig="",HNHSequence="",IMMSequence="",
                profiles=[],Species="",HNHstart=0,IMMstart=-1000000,HNHstop=-10000000,
                IMMstop=0, MATCH="FALSE", comments="", contig="", shortcontig="", contig_file=''):
      self.identifier = identifier
      self.coltype=coltype
      self.contiglength=contiglength
      self.database = database
      self.HNHcontig = HNHcontig
      self.IMMcontig=IMMcontig
      self.HNHsequence=HNHSequence
      self.HNHnucleotide=""
      self.IMMnucleotide=""
      self.IMMsequence=IMMSequence
      self.HNHstart=HNHstart
      self.IMMstart=IMMstart
      self.HNHstop=HNHstop
      self.IMMstop=IMMstop
      self.profiles=profiles
      self.Species=Species
      self.isolate=isolate
      self.MATCH=MATCH
      self.comments=comments
      self.VF=[]
      self.niche=""
      self.HNHORF=""
      self.IMMORF=""
      self.contig=contig
      self.shortcontig=shortcontig
      self.contig_end=0
      self.contig_start=0
      self.finder=""
      self.env_products=[]
      self.contig_file=''
      ISOLATE.isolateCount += 1

   def displayCount(self):
     print "Total sequences %d" % ISOLATE.isolateCount

   def displayISOLATE(self):
      print "Isolate : ", self.identifier,  ", startpos: ", self.HNHstart, ",endpos", self.HNHstop,"Species", self.Species

   def speciesWriter(self):
       string=self.isolate
       parts=string.split("_")
       string=str(parts[0]+"_"+parts[1])
       self.Species=string

   def contigCutter(self):
       if self.contig_start < self.contig_end:
           temp=self.contig_end
           self.contig_end=self.contig_start
           self.contig_start=temp
       if len(self.contig)>=30000:

           if self.contig_start-10000<0:
               alpha=0
           else:
               alpha=self.contig_start
           if self.contig_end +10000>len(self.contig):
               beta=0
           else:
               beta=self.contig_end

               self.shortcontig=self.contig[alpha:beta]
       else:
           self.shortcontig=self.contig
