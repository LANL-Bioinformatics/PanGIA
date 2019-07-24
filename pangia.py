#!/usr/bin/env python3
__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = "1.0.0-RC6.1"
__date__      = "2018/12/10"
__copyright__ = "BSD-3"

import argparse as ap
import textwrap as tw
import sys
import os
import time
import subprocess
import json
import gc
import glob
import gzip
from re import search, findall
from multiprocessing import Pool
from itertools import chain
import fileinput
import copy
import pandas as pd
import taxonomy as t
import pathogen as p
import score as s
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def parse_params(ver):
	class SmartFormatter(ap.HelpFormatter):
		def _split_lines(self, text, width):
			if text.startswith('R|'):
				return text[2:].splitlines()  
			# this is the RawTextHelpFormatter._split_lines
			return ap.HelpFormatter._split_lines(self, text, width)

	p = ap.ArgumentParser(prog='pangia.py', description="""PanGIA Bioinformatics %s""" % ver, formatter_class=SmartFormatter)

	eg = p.add_mutually_exclusive_group(required=True)

	eg.add_argument('-i', '--input', 
			metavar='[FASTQ]', nargs='+', type=str,
	  				help="Input one or multiple FASTQ file(s). Use space to separate multiple input files.")

	eg.add_argument('-s', '--sam', 
			metavar='[SAMFILE]', nargs=1, type=ap.FileType('r'),
					help="Specify the input SAM file. Use '-' for standard input.")

	p.add_argument('-d', '-db', '--database',
			metavar='[INDEX]', type=str, nargs='*',
					help="Name/path of BWA-MEM or Minimap2 index(es). [default: None]")

	p.add_argument('-dp', '--dbPath',
			metavar='[PATH]', type=str, default=None,
					help="""Path of databases. If this option isn't specified but a path is provided in "--database" option, this path of database will also be used in dbPath. 
					Otherwise, the program will search "database/" in program directory.
					[default: database/]""")

	default_val = 40
	p.add_argument('-asl','--alignSeedLength',
			metavar='<INT>', type=int, default=default_val,
			help="Minimum seed length uses in BWA-MEM [default: %s]"%default_val)

	default_val = 60
	p.add_argument('-ams','--alignMinScore',
			metavar='<INT>', type=int, default=default_val,
			help="Minimum alignment score (AS:i tag) for BWA-MEM [default: %s]"%default_val)

	default_val = '-h100 -B2'
	p.add_argument('-ao','--addOptions',
			metavar='<STR>', type=str, default=default_val,
			help="Additional options for BWA-MEM (no need to add -t) [default: '%s']"%default_val)

	p.add_argument('-rm', '--readmapper',
			choices=['bwa','minimap2'], type=str, default='bwa',
			help="Choose a read mapper from bwa or minimap2 [default: 'bwa']")

	p.add_argument('-se','--singleEnd', action="store_true",
					help="Input single-end reads or treat paired-end reads as single-end [default: False]")

	p.add_argument('-st','--scoreMethod',
			type=str, default='standalone',
					choices=['bg', 'standalone', 'combined'],
					help="""R|You can specify one of the following scoring method:\n"""
						 """"bg"         : compare mapping results with the background;\n"""
						 """"standalone" : score based on uniqueness;\n"""
						 """"combined"       : bg * standalone;\n"""
						 """[default: 'standalone']""" )

	p.add_argument('-sp','--scoreParam',
			type=str, default='0.5:0.99',
					help="""Normalization range of uniqueness score parameters [default: 0.5:0.99]""" )

	# p.add_argument('-sbgm','--scoreBgMethod',
	# 		type=str, default='overlapping_prop_raw',
	# 				choices=['overlapping_prop_raw', 'overlapping_prop_freq', 'KS_2samp_freq', 'KS_2samp_freq_removedZero'],
	# 				help="""R|You can specify one of the following scoring method:\n"""
	# 					 """"overlapping_prop_raw"       : overlapping proportion of raw depth distribution;\n"""
	# 					 """"overlapping_prop_freq"      : overlapping proportion of depth-frequency distribution;\n"""
	# 					 """"KS_2samp_freq"              : 2sample-KS-test depth-frequency distribution;\n"""
	# 					 """"KS_2samp_freq_removedZero"  : 2sample-KS-test depth-frequency distribution (depth=0 excluded);\n"""
	# 					 """[default: 'overlapping_prop_raw']""" )

	p.add_argument('-m','--mode',
			type=str, default='report',
					choices=['report', 'class', 'extract', 'lineage'],
					help="""R|You can specify one of the following output modes:\n"""
						 """"report"  : report a summary of profiling result;\n"""
						 """"class"   : output results of classified reads;\n"""
						 """"extract" : extract mapped reads;\n"""
						 """"lineage" : output abundance and lineage in a line;\n"""
						 """Note that only results/reads belongs to descendants of TAXID will be reported/extracted if option [--taxonomy TAXID] is specified. [default: 'report']""" )

	report_fields_choices = ['basic', 'r', 'rnb', 'rnr', 'ri', 'patho', 'score', 'ref', 'full', 'all']

	p.add_argument('-rf','--reportFields',
			type=str, default='all', nargs='+',
					choices=report_fields_choices,
					help="""R|You can specify following set of fields to display in the report:\n"""
						 """"basic" : essential fields that will display in the reports;\n"""
						 """"r"     : rank specific read count;\n"""
						 """"rnb"   : rank specific read count normalized by \n"""
						 """          both identity and # of ref (1*identity/num_refs);\n"""
						 """"rnr"   : rank specific read count normalized by \n"""
						 """          the number of references (1/num_refs);\n"""
						 """"ri"    : rank specific read identity\n"""
						 """          (mapped_length-nm)/read_length;\n"""
						 """"patho" : metadata of pathogen;\n"""
						 """"score" : detail score information;\n"""
						 """"ref"   : mapped reference(s) and their locations\n"""
						 """"full"  : display additional information\n"""
						 """"all"   : display all of above;\n"""
						 """[default: 'all']""")

	p.add_argument( '-da','--displayAll', action="store_true", required=False,
					help="Display all taxonomies including being filtered out [default: None]")

	p.add_argument( '-par','--procAltRefs', metavar='<INT>', type=int, default=30,
					help="Process the number of different references in alternative alignments [default: 30]")

	p.add_argument( '-xnm','--extraNM', metavar='<INT>', type=int, default=1,
					help="Process alternative alignments with extra number of mismatches than primary alignment [default: 1]")

	p.add_argument( '-x','--taxonomy', metavar='[TAXID]', type=str,
					help="Specify a NCBI taxonomy ID. The program  will only report/extract the taxonomy you specified.")

	p.add_argument( '-r','--relAbu', metavar='[FIELD]', type=str, default='DEPTH_COV',
					choices=['TOTAL_BP_MAPPED','READ_COUNT','PRI_READ_COUNT','DEPTH_COV','READ_COUNT_RNR','READ_COUNT_RSNB','RPKM'],
					help='The field will be used to calculate relative abundance. [default: DEPTH_COV]')

	p.add_argument( '-t','--threads', metavar='<INT>', type=int, default=1,
					help="Number of threads [default: 1]")

	p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
					help="Output directory [default: .]")

	p.add_argument( '-td','--tempdir', metavar='[DIR]', type=str,
					help="Default temporary directory [default: <OUTDIR>/<PREFIX>_tmp]")

	p.add_argument( '-kt','--keepTemp', action="store_true", required=False,
					help="Keep temporary directory after finishing the pipeline.")

	p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
					help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

	p.add_argument( '-ps','--pathoScoreOnly', action="store_true", required=False,
					help="Only calculate score for pathogen under '--scoreMethod bg'")

	p.add_argument( '-sb','--saveBg', action="store_true", required=False,
					help="Save current readmapping result in JSON to <PREFIX>.json")

	p.add_argument( '-lb','--loadBg', metavar='<FILE>', nargs='*', type=str, required=False,
					help="Load one or more background JSON gzip file(s) [default: None")

	default_val = 0
	p.add_argument( '-ms','--minScore', metavar='<FLOAT>', type=float, default=default_val,
					help="Minimum score to be considered valid [default: %s]"%default_val)

	default_val = 10
	p.add_argument( '-mr','--minReads', metavar='<INT>', type=int, default=default_val,
					help="Minimum number of reads to be considered valid [default: %s]"%default_val)

	default_val = 2.5
	p.add_argument( '-mb','--minRsnb', metavar='<INT>', type=int, default=default_val,
					help="Minimum number of reads to be considered valid [default: %s]"%default_val)

	default_val = 200
	p.add_argument( '-ml','--minLen', metavar='<INT>', type=int, default=default_val,
					help="Minimum linear length to be considered valid [default: %s]"%default_val)

	default_val = 0.004
	p.add_argument( '-mc','--minCov', metavar='<FLOAT>', type=float, default=default_val,
					help="Minimum linear coverage to be considered a valid strain [default: %s]"%default_val)

	default_val = 0.01
	p.add_argument( '-md','--minDc', metavar='<FLOAT>', type=float, default=default_val,
					help="Minimum depth of coverage to be considered a valid strain [default: %s]"%default_val)

	default_val = 0.0009
	p.add_argument( '-mrd','--minRsdcnr', metavar='<FLOAT>', type=float, default=default_val,
					help="Minimum rank specific depth of coverage normalized by the number of mapped references to be considered a valid strain [default: %s]"%default_val)

	# p.add_argument( '-np','--nanopore', action="store_true",
	# 				help="Input reads is nanopore data. This option is equivalent to use [-oa='-h 150 -x ont2d' -ms 0 -mr 1 -mb 3 -ml 50 -asl 24 -ams 70]. [default: FALSE]")
	p.add_argument( '-np','--nanopore', action="store_true",
					help="Input nanopore reads and use minimap2 to align reads. [default: FALSE]")

	p.add_argument( '-pd','--pathogenDiscovery', action="store_true",
					help="Adjust options for pathogen discovery. This option is equivalent to use [-ms 0 -mr 3 -mb 1 -ml 50 -asl 24 -ams 50 -mc 0 -md 0 -mrd 0]. [default: FALSE]")

	p.add_argument( '-if','--ignoreFlag', metavar='<STR>', type=str, default=None,
					help="Ignore reads that mapped to the references that have the flag(s) [default: None]")

	p.add_argument( '-nc','--noCutoff', action="store_true",
					help="Remove all cutoffs. This option is equivalent to use [-ms 0 -mr 0 -mb 0 -ml 0 -mc 0 -md 0 -mrd 0].")

	p.add_argument( '-c','--stdout', action="store_true",
					help="Write on standard output.")

	p.add_argument( '--silent', action="store_true",
					help="Disable all messages.")

	p.add_argument( '--verbose', action="store_true",
					help="Provide verbose running messages and keep all temporary files.")

	p.add_argument( '--debug', action="store_true",
					help="Provide debugging messages and keep all temporary files.")

	eg.add_argument( '--version', action="store_true", default=False,
					help="Print version number.")

	args_parsed = p.parse_args()

	"""
	Checking options
	"""
	if args_parsed.version:
		print( ver )
		os._exit(0)

	if args_parsed.input and not args_parsed.database:
		p.error( '--database option is missing.' )

	if args_parsed.input and args_parsed.sam:
		p.error( '--input and --same are incompatible options.' )

	if not args_parsed.dbPath:
		# set dbPath to the path of the first input database if the input path is absolute path
		if args_parsed.database and "/" in args_parsed.database[0]:
			db_dir = search( '^(.*?)[^\/]+$', args_parsed.database[0] )
			args_parsed.dbPath = db_dir.group(1)
		else:
			bin_dir = os.path.dirname(os.path.realpath(__file__))
			args_parsed.dbPath = bin_dir + "/database"

	if not os.path.isfile( args_parsed.dbPath + "/taxonomy.tsv" ) and not os.path.isfile( args_parsed.dbPath + "/names.dmp" ):
		p.error( "NCBI dmp files or taxonomy.tsv not found at %s. Please check --dbPath option."%args_parsed.dbPath )
	if not os.path.isfile( args_parsed.dbPath + "/taxonomy.custom.tsv" ):
		p.error( "taxonomy.custom.tsv not found. Please check --dbPath option." )
	if not os.path.isfile( args_parsed.dbPath + "/pathogen.tsv" ):
		p.error( "pathogen.tsv not found. Please check --dbPath option." )
	if not os.path.isfile( args_parsed.dbPath + "/uniqueness.tsv" ):
		p.error( "uniqueness.tsv not found. Please check --dbPath option." )

	if args_parsed.scoreMethod != "standalone":
		if not args_parsed.loadBg:
			p.error( 'Missing background info: please specify with --loadBg option.' )

	if "all" in args_parsed.reportFields:
		args_parsed.reportFields = report_fields_choices

	if not args_parsed.prefix:
		if args_parsed.input:
			name = search('([^\/\.]+)\..*$', args_parsed.input[0] )
			args_parsed.prefix = name.group(1)
		elif args_parsed.sam:
			name = search('([^\/]+).\w+.\w+$', args_parsed.sam[0].name )
			args_parsed.prefix = name.group(1)
		else:
			args_parsed.prefix = "pangia"

	if not args_parsed.tempdir:
		args_parsed.tempdir = args_parsed.outdir+"/"+args_parsed.prefix+"_tmp"

	if not args_parsed.singleEnd:
		args_parsed.singleEnd = "auto"

	if args_parsed.debug:
		args_parsed.verbose = True

	if args_parsed.verbose or args_parsed.debug:
		args_parsed.keepTemp = True

	if args_parsed.readmapper == "minimap2":
		args_parsed.addOptions = "-A1 -B2 -k %s -m %s -x sr -p 1 -N 30"%(args_parsed.alignSeedLength, args_parsed.alignMinScore)

	if args_parsed.nanopore:
		args_parsed.readmapper = "minimap2"
		args_parsed.addOptions = "-x map-ont"
		args_parsed.minScore = 0
		args_parsed.minRsnb = 1
		args_parsed.minReads = 3
		args_parsed.minLen = 50

	if args_parsed.pathogenDiscovery:
		args_parsed.minScore = 0
		args_parsed.minRsnb = 1
		args_parsed.minReads = 3
		args_parsed.minLen = 50
		args_parsed.alignSeedLength = 24
		args_parsed.alignMinScore = 50
		args_parsed.minCov = 0
		args_parsed.minDc = 0
		args_parsed.minRsdcnr = 0

	if args_parsed.noCutoff:
		args_parsed.minScore = 0
		args_parsed.minRsnb = 0
		args_parsed.minReads = 0
		args_parsed.minLen = 0
		args_parsed.minCov = 0
		args_parsed.minDc = 0
		args_parsed.minRsdcnr = 0

	# glob database path
	if args_parsed.database:
		new_db=[]
		for db in args_parsed.database:
			if args_parsed.readmapper == "bwa":
				if args_parsed.dbPath and not "/" in db and not os.path.isfile( db+".ann" ):
					db = args_parsed.dbPath+"/"+db
				if not os.path.isfile( db + ".ann" ):
					p.error( 'Incorrect BWA index: missing %s.ann.' % db )
			elif args_parsed.readmapper == "minimap2":
				# the input database for minimap2 can be a raw fasta file or pre-build index (*.mmi).
				if args_parsed.dbPath and not "/" in db:
					db = args_parsed.dbPath+"/"+db
				if os.path.isfile( db + ".mmi" ):
					db += ".mmi"
			
			new_db.append(db)
		args_parsed.database = new_db

	return args_parsed

def dependency_check(cmd):
	""" Check whether `cmd` is on PATH and executable. """
	# from whichcraft import which
	from shutil import which
	cmd_path = which(cmd)
	return cmd_path if cmd_path is not None else False

def isDescendant( taxid, taxid_ant ):
	fullLineage = t.taxid2fullLineage( taxid )
	if "|%s|" % taxid_ant in fullLineage:
		return True
	else:
		return False

def lineageLCR(taxids):
	""" lineageLCR
	
	Find late common rank (LCR) from two lineage

	Arguments:
	  arg1 <LIST>: taxids

	Returns:
	  level <STR>: LCR name of current and quering lineage
	  name <STR> : LCR level of current and quering lineage
	"""

	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']

	merged_dict = t._autoVivification()
	for tid in taxids:
		lng = t.taxid2lineageDICT(tid, 1, 1)
		for r in lng:
			ttid = lng[r]['taxid']
			if ttid in merged_dict[r]:
				merged_dict[r][ttid] += 1
			else:	
				merged_dict[r][ttid] = 1

	for r in ranks:
		if len(merged_dict[r]) == 1:
			for ttid in merged_dict[r]:
				# skip if no tid in this rank
				if ttid==0:
					continue
				tname = t.taxid2name(ttid)
				return r, tname, merged_dict

	return "root", "root", merged_dict

def flagInHeader( flag, header ):
	if flag:
		h = header.split('|')
		for f in flag:
			if f in h[3]:
				return True

	return False

def parse(line, ignore_flag, procAltRefs, extraNM):
	"""Parse SAM file

	Example input:
	@SQ     SN:ref1 LN:40
	@SQ     SN:ref2 LN:40
	read1	0	ref2	11	0	10M	*	0	0	GGGGGGGGGG	*	NM:i:0	MD:Z:10	AS:i:10	XS:i:10	XA:Z:ref3,+21,10M,0;ref1,+1,10M,0;ref4,+31,10M,0;
	read2	16	ref4	31	0	10M	*	0	0	GGGGGGGGGG	*	NM:i:0	MD:Z:10	AS:i:10	XS:i:10	XA:Z:ref1,-1,10M,0;ref3,-21,10M,0;ref2,-11,10M,0;	
	XA tag contains suboptional alignments in (chr,pos,CIGAR,NM;) format

	Arguments:
	  arg1 <STR>: alignment result of a read in SAM format

	Returns:
	  value <DICT>:
		-> 'name' <STR>: read name 
		-> 'seq' <STR> : read seqeunce
		-> 'ref' <DICT>: target referecnes 
		  -> [ref name] <DICT>: referecne header,
							   for example: "NC00001|766|23762|Bc"
			-> 'start' <INT>  : mapping start position
			-> 'end'   <INT>  : mapping stop position <INT>
			-> 'nm'    <INT>  : number of mismatch <INT>
			-> 'idnt'  <FLOAT>: identity --> (mapped_len-mismatch_len)/read_len
			-> 'rlen'  <INT>  : reference length <INT>
			-> 'std'   <INT>  : strand + or -
			-> 'taxid' <STR>  : taxonomy ID <STR>
	"""

	#prepare output file
	read_mr = t._autoVivification() #read mapping result

	#split SAM result
	temp   = line.split('\t')
	read_mr['name'] = temp[0]
	samflag = int(temp[1])
	read_mr['raw'] = line

	#ignore read if flag found
	if flagInHeader( ignore_flag, temp[2]):
		return {'name':temp[0]}

	#read reverse strand
	if samflag & 16:
		read_mr['seq'] = seqReverseComplement(temp[9])
		read_mr['qua'] = temp[10][::-1]
	else:
		read_mr['seq'] = temp[9]
		read_mr['qua'] = temp[10]

	#first or second mate
	if samflag & 64:
		read_mr['pe'] = 1
	elif samflag & 128:
		read_mr['pe'] = 2	

	#alignment score (primary)
	s_as = search('AS:i:(\d+)', line)
	#DP score of the max scoring segment in the alignment (minimap2)
	s_ms = search('ms:i:(\d+)', line)
	try:
		read_mr['score'] = int(s_as.group(1))
	except:
		try:
			read_mr['score'] = int(s_ms.group(1))
		except:
			read_mr['score'] = 0

	# remove the alignment if < alignMinScore
	if read_mr['score'] < argvs.alignMinScore:
		return {'name':temp[0]}
	
	#alignment score (alternative)
	s_xs = search('XS:i:(\d+)', line)
	#minimap2
	try:
		read_mr['xs'] = int(s_xs.group(1))
	except:
		read_mr['xs'] = 0

	#skip unmapped reads
	if int(temp[1]) & 4: return {'name':temp[0]}

	#primary alignment
	read_len       = len(temp[9])
	s_mismatch_len = search('NM:i:(\d+)', line)
	mismatch_len   = int(s_mismatch_len.group(1))
	clipped_len    = sum([ int(s) for s in findall('(\d+)S', temp[5]) ])
	unclipped_len  = read_len-clipped_len
	mapped_len     = sum([ int(s) for s in findall('(\d+)M', temp[5]) ])
	
	start          = int(temp[3])
	end            = start + unclipped_len - 1

	# skip alignment if "unknown" taxonomy id found  
	if not '|unknown|' in temp[2]:
		read_mr['ref'][temp[2]]['start'] = start
		read_mr['ref'][temp[2]]['end']   = end
		read_mr['ref'][temp[2]]['nm']    = mismatch_len
		read_mr['ref'][temp[2]]['idnt']  = (mapped_len-mismatch_len)/read_len
		#primary alignment flag
		read_mr['ref'][temp[2]]['p']     = True
		#strand
		read_mr['ref'][temp[2]]['strand'] = '-' if int(temp[1]) & 16 else '+'

	#process alternative alignments that has >80% of the best score (BWA-MEM built-in)
	if temp[-1].startswith('XA'):
		xa_tag = search('XA:Z:(\S+);', temp[-1])
		alns = xa_tag.group(1).split(';')

		# iterating over XAs
		# Only the first <procAltRefs> number of references will be processed due to performance (sorting a large number of XAs is very slow).
		for aln in alns[0:procAltRefs]:
			attrs = aln.split(',')
			ref_n = attrs[0]
			ref_nm = int(attrs[3])
			mapped_len = sum([ int(s) for s in findall('(\d+)M', attrs[2]) ])

			# Allow 1 more mismatch compared to the best score
			# BWA mem will report score > 80% of the best score, therefore,
			# if the NM of secondery alignment is more than the NM of the primary alignment,
			# the alignment will be filtered out (not the best alignment).
			if ref_nm > mismatch_len+extraNM: continue

			# skip reference with unknown taxonomy id 
			if '|unknown|' in ref_n: continue

			#only record one alignment in a replicon
			if ref_n in read_mr['ref']: continue

			#ignore read if flag found
			if flagInHeader( ignore_flag, ref_n): return {'name':temp[0]}
				
			clipped_len   = sum([ int(s) for s in findall('(\d+)S', attrs[2]) ])
			unclipped_len = read_len-clipped_len
			start         = int(attrs[1].lstrip('+-'))
			end           = start+unclipped_len-1

			read_mr['ref'][ref_n]['XA']     = aln
			read_mr['ref'][ref_n]['start']  = start
			read_mr['ref'][ref_n]['end']    = end
			read_mr['ref'][ref_n]['nm']     = ref_nm
			read_mr['ref'][ref_n]['idnt']   = (unclipped_len-ref_nm)/read_len
			read_mr['ref'][ref_n]['strand'] = '-' if '-' in attrs[1] else '+'
	else:
		# If score of secondary hit >= primary hit but no XA tags, discard this alignment.
		# Because this read hit >30 locations (by default), it likely a repeat or some low complexity region 
		if read_mr['score'] <= read_mr['xs']:
			return {'name':temp[0]}

	#parse ref name to get length and taxid
	for ref_n in read_mr['ref']:
		rs = ref_n.split('|')
		read_mr['ref'][ref_n]['rlen'] = int(rs[1])
		read_mr['ref'][ref_n]['taxid'] = rs[2]

	return read_mr

def worker(filename, chunkStart, chunkSize):
	"""worker(lines)

	Take alignment results to a big dict indexed by reference name. 

	Arguments:
	  arg1 <STR>: alignment result of a read in SAM format
	
	Return:
	  A dict mapping keys (ref header) to corresponding data
	  result <DICT>: read mapping results
		-> ref_n <DICT>: reference header
		  -> "ML" <LIST>: list of range in list [start,end]
		  -> "MB" <INT> : # of mapped bases (mapping region only -> match + mismatch)
		  -> "TB" <INT> : total length of mapped reads (match + mismatch) <INT>
		  -> "MR" <INT> : # of mapped reads <INT>
		  -> "NM" <INT> : # of mismatches <INT>
		  -> "LL" <INT> : linear length
		  -> "[lcr]" <STR>: #LCR, can be strain, species...up to superkingdom
		   			-> "R"     <INT>  : raw read count
					-> "M"     <INT>  : mapped length (match + mismatch)
					-> "Mnr"   <INT>  : mapped length / num_refs
					-> "T"     <INT>  : total read length
					-> "N"     <INT>  : total mismatch bases
					-> "Rnb"   <FLOAT>: normalize read count by both identity and # of ref (Ri/num_refs)
					-> "Rnr"   <FLOAT>: normalize read count by the number of references (1/num_refs)
					-> "Ri"    <FLOAT>: identity (mapped_length-nm)/read_length
					-> "R_IDT" <LIST> : list of the identity of each read
					-> "R_MAT" <LIST> : list of the # of matches of each read
					-> "R_SCR" <LIST> : list of alignment scores
	"""
	
	res = t._autoVivification()
	read_buffer = [] #keep 2 reads
	out_buffer = {}

	#benchmark start time
	stime = time.time()
	proced_num = 0
	def parsingBenchmark(stime, cnt):
		if (time.time()-stime)%10 < 0.02:
			nowtime = time.time()
			print( "%s\t%.2f\t%.2f" % (cnt, float(nowtime-stime), float(cnt/(nowtime-stime))) )

	def remove_unpaired_refs(read_buffer, r2):
		"""
		Delete references with no proper paired-read
		"""
		r1 = read_buffer[0]

		# return empty buffer if no ref found in the alignment
		if not 'ref' in r1 or not 'ref' in r2: return []

		if r1['name'] != r2['name']:
			return []
		else:
			#purage r1 refs
			all_ref_n = list(r1['ref'].keys())
			
			for ref_n in all_ref_n:
				if not ref_n in r2['ref']:
					del r1['ref'][ref_n]
				elif r1['ref'][ref_n]['strand'] == r2['ref'][ref_n]['strand']:
					del r1['ref'][ref_n]
					del r2['ref'][ref_n]
			#purage r2 refs
			all_ref_n = list(r2['ref'].keys())
			for ref_n in all_ref_n:
				if not ref_n in r1['ref']:
					del r2['ref'][ref_n]

		if len(r2['ref']) > 0:
			read_buffer.append(r2)
		else:
			return []

		return read_buffer

	def prepAlnForOutput(aln, out_buffer):
		"""
		This function expends XA tags in a SAM alignment to multiple SAM entries. The expanded
		SAM entries will append to out_buffer using reference name as keys.
		"""

		for ref_n in aln['ref']:
			# init output sam buffer for a reference
			if not ref_n in out_buffer: out_buffer[ref_n] = ''

			t = aln['raw'].split('\t')

			# print SAM out if it's primary reference
			if 'p' in aln['ref'][ref_n]:
				out_buffer[ref_n] += '\t'.join(t[:9]) + '\t*\t*\t' + '\t'.join(t[11:-1]) + '\n'
				continue

			a = aln['ref'][ref_n]['XA'].split(',')
			mchr = '=' if t[6] == a[0] else t[6] # FIXME: TLEN/ISIZE is not calculated!
			seq="*"
			phred="*"
			out_buffer[ref_n] += "\t".join( [ t[0], str((int(t[1])&0x6e9)|(0x10 if int(a[1])<0 else 0)), a[0], str(abs(int(a[1]))), "0", a[2], t[6], t[7], "0", seq, phred, "NM:i:"+a[3] ]) + "\n"

		return out_buffer

	# processing alignments in SAM format
	f = open( filename )
	f.seek(chunkStart)
	lines = f.read(chunkSize).splitlines()
	for line in lines:
		# ignore SAM header
		if line.startswith('@'): continue
		# parse alignment output in SAM format
		aln = parse(line, argvs.ignoreFlag, argvs.procAltRefs, argvs.extraNM )

		# paired-end or single-end reads
		if argvs.singleEnd == "auto":
			tmp = line.split("\t")
			if int(tmp[1]) & 1:
				argvs.singleEnd = False
			else:
				argvs.singleEnd = True

		if not argvs.singleEnd:
			# if paired end buffer is empty, push aln to the buffer
			if not read_buffer:
				read_buffer.append(aln)
				continue
			# if previous aln is in the buffer, 
			else:
				# if the alignment in the buffer doesn't belong to current read,
				# clean up buffer and push current aln in the buffer
				if read_buffer[0]['name'] != aln['name']:
					read_buffer = []
					read_buffer.append(aln)
					continue
				else:		
					read_buffer = remove_unpaired_refs(read_buffer, aln)
		else:
			read_buffer = []
			read_buffer.append(aln)

		for aln in read_buffer:
			# skip unmapped and ignored (post-processed) reads
			if not 'ref' in aln: continue

			#find LCR of 
			taxids = []
			for ref_n in aln['ref']:
				taxids.append(aln['ref'][ref_n]['taxid'])

			num_refs = len(aln['ref'])
			lcr_lvl, lcr_name, lcr_info = lineageLCR(taxids)

			# output expanded SAM to files if its specific enough
			if not lcr_lvl in ('class','phylum','superkingdom'):
				out_buffer = prepAlnForOutput(aln, out_buffer)

			for ref_n in aln['ref']:
				start  = aln['ref'][ref_n]['start']
				end    = aln['ref'][ref_n]['end']
				r      = [start, end]
				mapped_len = end - start + 1
				rlen   = len(aln['seq'])
				nm     = aln['ref'][ref_n]['nm']
				pri_read = True if 'p' in aln['ref'][ref_n] else False

				if ref_n in res:
					res[ref_n]["MB"] += mapped_len
					res[ref_n]["TB"] += rlen
					res[ref_n]["MR"] += 1
					res[ref_n]["NM"] += nm
					if pri_read: res[ref_n]["PR"] += 1
				else:
					res[ref_n]["MB"] = mapped_len
					res[ref_n]["TB"] = rlen
					res[ref_n]["MR"] = 1
					res[ref_n]["NM"] = nm
					res[ref_n]["PR"] = 1 if pri_read else 0

				if not lcr_lvl in res[ref_n]:
					res[ref_n][lcr_lvl]["R"]     = 0 
					res[ref_n][lcr_lvl]["M"]     = 0
					res[ref_n][lcr_lvl]["Mnr"]   = 0 
					res[ref_n][lcr_lvl]["T"]     = 0 
					res[ref_n][lcr_lvl]["N"]     = 0 
					res[ref_n][lcr_lvl]["Rnb"]   = 0 
					res[ref_n][lcr_lvl]["Rnr"]   = 0 

				res[ref_n][lcr_lvl]["R"]     += 1 #raw read count at rank
				res[ref_n][lcr_lvl]["M"]     += mapped_len
				res[ref_n][lcr_lvl]["Mnr"]   += mapped_len/num_refs
				res[ref_n][lcr_lvl]["T"]     += rlen
				res[ref_n][lcr_lvl]["N"]     += nm
				res[ref_n][lcr_lvl]["Rnb"]   += (mapped_len-nm)/rlen/num_refs #normalize read count by references and identity
				res[ref_n][lcr_lvl]["Rnr"]   += 1/num_refs #normalize read count by the number of references

	# create a directory to save sam files for this subprocess
	if not os.path.exists( "%s/sam_%s"%(argvs.tempdir, os.getpid()) ):
		os.makedirs( "%s/sam_%s"%(argvs.tempdir, os.getpid()) )
	# close output file handelers
	for ref_n in out_buffer:
		filename = "%s/sam_%s/%s.sam"%(argvs.tempdir, os.getpid(), ref_n)
		with open( filename, "a") as f:
			f.write( out_buffer[ref_n] )
			f.close()

	return res

def linear_cov_calculation( tempdir, threads ):
	"""
	linear_cov_calculation
	"""

	lc={}
	script_content = """#!/bin/bash
set -e
OUTDIR=$1
THREAD=$2
cd $OUTDIR
find sam_* -type f -name '*.sam' | rev | cut -d'/' -f1 | sort | uniq | rev > sam_uniq.list
mkdir -p merged_sam
parallel -j$THREAD --halt 2 'printf "@SQ\\tSN:"%s"\\tLN:"%s"\\n" `printf {/.}` `printf {/.} | cut -d"|" -f2` > merged_sam/{};\
  cat sam_*/{} >> merged_sam/{};\
  samtools view -bhS merged_sam/{} > merged_sam/{/.}.unsorted.bam;\
  samtools sort merged_sam/{/.}.unsorted.bam > merged_sam/{/.}.sorted.bam;\
  samtools index merged_sam/{/.}.sorted.bam;\
  samtools depth -d1000000 merged_sam/{/.}.sorted.bam > merged_sam/{/.}.sorted.depth;\
  printf {/.}"\\t";\
  cat merged_sam/{/.}.sorted.depth | wc -l' :::: sam_uniq.list > ref_linear_cov.txt
set +e
"""
	with open('%s/ref_linear_cov.sh'%tempdir, 'w') as s:
		s.write(script_content)
		s.close()

	cmd = "set -o pipefail; set -x; bash %s/ref_linear_cov.sh %s %s"%(tempdir, tempdir, threads)

	if argvs.verbose: print_message( "[INFO] CMD: %s"%cmd, argvs.silent, begin_t, logfile )

	proc = subprocess.Popen( cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	outs, errs = proc.communicate()
	exitcode = proc.poll()

	if exitcode == 0:
		with open('%s/ref_linear_cov.txt'%tempdir) as c:
			for line in c:
				(ref_n, ll) = line.rstrip().split('\t') #get reference name and linear length
				lc[ref_n] = int(ll) 
			c.close()
	else:
		print_message( "[%s] ERROR: error occurred while calculating linear coverage (exit: %s, message: %s).\n" % (timeSpend(begin_t), exitcode, errs), argvs.silent, begin_t, logfile, True )

	return lc

def result_merging_worker( result_list ):
	result = t._autoVivification()
	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom','root']

	# merge all replicon results processed by worker function in parallel
	for res in result_list:
		for k in res:
			if k in result:
				result[k]["MB"] += res[k]["MB"]
				result[k]["TB"] += res[k]["TB"]
				result[k]["MR"] += res[k]["MR"]
				result[k]["PR"] += res[k]["PR"]
				result[k]["NM"] += res[k]["NM"]
				for lcr_lvl in ranks:
					if lcr_lvl in res[k]:
						if not lcr_lvl in result[k]:
							result[k][lcr_lvl]["R"]      = 0
							result[k][lcr_lvl]["M"]      = 0
							result[k][lcr_lvl]["Mnr"]    = 0
							result[k][lcr_lvl]["T"]      = 0
							result[k][lcr_lvl]["N"]      = 0
							result[k][lcr_lvl]["Rnb"]    = 0
							result[k][lcr_lvl]["Rnr"]    = 0

						result[k][lcr_lvl]["R"]      += res[k][lcr_lvl]["R"]
						result[k][lcr_lvl]["M"]      += res[k][lcr_lvl]["M"]
						result[k][lcr_lvl]["Mnr"]    += res[k][lcr_lvl]["Mnr"]
						result[k][lcr_lvl]["T"]      += res[k][lcr_lvl]["T"]
						result[k][lcr_lvl]["N"]      += res[k][lcr_lvl]["N"]
						result[k][lcr_lvl]["Rnb"]    += res[k][lcr_lvl]["Rnb"]
						result[k][lcr_lvl]["Rnr"]    += res[k][lcr_lvl]["Rnr"]
			else:
				result[k].update(res[k])
	return result

def chunkify(fname, size=1*1024*1024):
	fileEnd = os.path.getsize(fname)
	with open(fname, "rb") as f:
		chunkEnd = f.tell()
		while True:
			chunkStart = chunkEnd
			f.seek(size, 1)
			f.readline()
			# put paired-end reads in the same chunck
			line = f.readline().decode('ascii')
			if chunkEnd <= fileEnd and line:
				tmp = line.split('\t')
				if int(tmp[1]) & 1 and int(tmp[1]) & 64:
					f.readline()
			# current position
			chunkEnd = f.tell()
			yield chunkStart, chunkEnd - chunkStart
			if chunkEnd > fileEnd:
				break

def processSAMfile( sam_fn, numthreads, numlines ):
	"""processSAMfile

	Read read mapping results in SAM foramt and return a dict.

	Arguments:
	  sam_fn <STR>    : filename of SAM file
	  numthreads <INT>: number of threads
	  numlines <INT>  : split every # of lines
	
	Returns:
	  result <DICT>: read mapping results
		  -> ref_n <DICT>: reference header as key
			  -> "ML" <LIST>: list of range in list [start,end]
			  -> "MB" <INT> : # of mapped bases (mapping region only -> match + mismatch)
			  -> "TB" <INT> : total length of mapped reads (match + mismatch) <INT>
			  -> "MR" <INT> : # of mapped reads <INT>
			  -> "NM" <INT> : # of mismatches <INT>
			  -> "LL" <INT> : linear length
			  -> rank <DICT>: # LCR rank
		   			-> "R"     <INT>  : raw read count
					-> "M"     <INT>  : mapped length (match + mismatch)
					-> "T"     <INT>  : total read length
					-> "N"     <INT>  : total mismatch bases
					-> "Rnb"   <FLOAT>: normalize read count by both identity and # of ref (1*identity/num_refs)
					-> "Rnr"   <FLOAT>: normalize read count by the number of references (1/num_refs)
					-> "Rni"   <FLOAT>: identity (mapped_length-nm)/read_length
	"""

	result = t._autoVivification()
	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom','root']

	#clean memory
	gc.collect()

	print_message( "Parsing SAM files with %s subprocesses..."%numthreads, argvs.silent, begin_t, logfile )

	if numthreads > 1:
		pool = Pool(processes=numthreads)
		jobs = []
		results = []

		for chunkStart,chunkSize in chunkify(sam_fn):
			jobs.append( pool.apply_async(worker, (sam_fn,chunkStart,chunkSize)) )

		#wait for all jobs to finish
		tol_jobs = len(jobs)
		cnt=0
		for job in jobs:
			results.append( job.get() )
			cnt+=1
			if argvs.verbose: print_message( "[INFO] Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

		#clean up
		pool.close()

		print_message( "Merging results...", argvs.silent, begin_t, logfile )
		result = result_merging_worker(results)
	else:
		result = worker(sam_fn, 0, os.path.getsize(sam_fn))

	print_message( "Done.", argvs.silent, begin_t, logfile )

	# Add linear length to results
	print_message( "Calculating linear length...", argvs.silent, begin_t, logfile )
	lc = linear_cov_calculation( argvs.tempdir, argvs.threads )

	tol_num_mapped_reads = 0
	for k in result:
		if k in lc:
			result[k]["LL"] = lc[k]
		else:
			result[k]["LL"] = 0
		tol_num_mapped_reads += result[k]["PR"]

	return result, tol_num_mapped_reads

def processSAMfileReadClass( f, o, taxid ):
	"""processSAMfileReadClass

	This function is used to extract reads that belongs to specified taxonomy ID.

	Arguments:
		f <FILE>: SAM file
		o <FILE>: output FASTQ file
		taxid <STR>: user specified taxonomy ID

	Returns:
		None
	"""

	for line in f:
		# parse alignment output in SAM format
		sam = parse(line, argvs.ignoreFlag, argvs.procAltRefs, argvs.extraNM )

		# continue if unmapped
		if not 'ref' in sam: continue

		# check if any ref is belong to provided taxid
		if taxid:
			in_descendant = False
			for ref_n in sam['ref']:
				if isDescendant( sam['ref'][ref_n]['taxid'], taxid ):
					in_descendant = True
					break
		else:
			in_descendant = True

		if not in_descendant: continue 

		#find LCR of alignment and location
		location = ""
		taxids = []
		for ref_n in sam['ref']:
			taxids.append(sam['ref'][ref_n]['taxid'])
			location += "%s:%s..%s;" % (ref_n, sam['ref'][ref_n]['start'], sam['ref'][ref_n]['end'])

		lcr_lvl, lcr_name, lcr_info = lineageLCR(taxids)

		o.write( "%s\t[%s|%s]\t%s\n" % (
			sam['name'],
			lcr_lvl,
			lcr_name,
			location
		))

def processSAMfileReadExtract( sam_fn, o, taxid, numthreads ):
	pool = Pool(processes=numthreads)
	jobs = []

	for chunkStart,chunkSize in chunkify(sam_fn):
		jobs.append( pool.apply_async(ReadExtractWorker, (sam_fn,chunkStart,chunkSize,taxid)) )

	#wait for all jobs to finish
	for job in jobs:
		outread = job.get()
		o.write(outread)
		o.flush()

	#clean up
	pool.close()

def ReadExtractWorker( filename, chunkStart, chunkSize, taxid ):
	# output
	readstr=""
	# processing alignments in SAM format
	f = open( filename )
	f.seek(chunkStart)
	lines = f.read(chunkSize).splitlines()
	for line in lines:
		# parse alignment output in SAM format
		sam = parse(line, argvs.ignoreFlag, argvs.procAltRefs, argvs.extraNM )

		# continue if unmapped
		if not 'ref' in sam: continue

		# check if any ref is belong to provided taxid
		if taxid:
			in_descendant = False
			for ref_n in sam['ref']:
				if isDescendant( sam['ref'][ref_n]['taxid'], taxid ):
					in_descendant = True
					break
		else:
			in_descendant = True

		if not in_descendant: continue

		#find LCR of alignment and location
		location = ""
		taxids = []
		for ref_n in sam['ref']:
			taxids.append(sam['ref'][ref_n]['taxid'])
			location += "%s:%s..%s;" % (ref_n, sam['ref'][ref_n]['start'], sam['ref'][ref_n]['end'])

		lcr_lvl, lcr_name, lcr_info = lineageLCR(taxids)

		readstr += "@%s  [%s|%s]  %s\n%s\n+\n%s\n" % (
			sam['name'],
			lcr_lvl,
			lcr_name,
			location,
			sam['seq'],
			sam['qua']
		)
	return readstr

def seqReverseComplement( seq ):
	"""Reverse complement nucleotide sequence

	Arguments:
	  arg1 <STR>: nucleotide sequence

	Return:
	  <STR>: reverse-complement input requence
	"""
	for base in seq:
		if base not in 'ACGTURYSWKMBDHVNacgturyswkmbdhvn':
			print("[ERROR] NOT a DNA sequence")
			return None
	seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
	seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
	return "".join([seq_dict[base] for base in reversed(seq)])

def taxonomyRollUp(r, patho_meta, tol_num_mapped_reads, mb, mr, ml, mc, md):
	"""taxonomyRollUp
	
	Take parsed SAM output and rollup to superkingdoms

	Arguments:
	arg1:
	  r <DICT>: read mapping results
		  -> ref_n <DICT>: reference header as key
			  -> "ML" <LIST>: list of range in list [start,end]
			  -> "MB" <INT> : # of mapped bases (mapping region only -> match + mismatch)
			  -> "TB" <INT> : total length of mapped reads (match + mismatch) <INT>
			  -> "MR" <INT> : # of mapped reads <INT>
			  -> "PR" <INT> : # of primary mapped reads <INT>
			  -> "NM" <INT> : # of mismatches <INT>
			  -> "LL" <INT> : linear length
			  -> rank <DICT>: # LCR rank
		   			-> "R"     <INT>  : raw read count
					-> "M"     <INT>  : mapped length (match + mismatch)
					-> "T"     <INT>  : total read length
					-> "N"     <INT>  : total mismatch bases
					-> "Rnb"   <FLOAT>: normalize read count by both identity and # of ref (1*identity/num_refs)
					-> "Rnr"   <FLOAT>: normalize read count by the number of references (1/num_refs)
					-> "Ri"    <FLOAT>: identity (mapped_length-nm)/read_length

	Returns:
	  res_rollup <DICT>: read mapping results roll up to superkingdom
		  -> taxid <DICT>: taxid as keys
			  -> "ML"   <LIST>: list mapped reference(s) and their locations
			  -> "MB"   <INT> : # of mapped bases (mapping region only -> match + mismatch)
			  -> "TB"   <INT> : total length of mapped reads (match + mismatch)
			  -> "MR"   <INT> : # of mapped reads
			  -> "PR"   <INT> : # of primary mapped reads
			  -> "DC"   <FLOAT> : depth of coverage
			  -> "DCnr" <FLOAT> : depth of coverage normalized by the number of references
			  -> "rsDCnr" <FLOAT> : rank specific depth of coverage normalized by the number of references
			  -> "MRNr" <INT> : Total normalize read count by the number of references
			  -> "rsMRNb" <INT> : Total rank specific normalized reads
			  -> "rsMB" <INT> : Total rank specific mapped bases
			  -> "RPKM" <INT> : RPKM
			  -> "NM"   <INT> : # of mismatches
			  -> "LL"   <INT> : linear length (bp)
			  -> "LC"   <INT> : best linear coverage (prop)
			  -> "P"    <STR> : Pathogen taxonomy id; 0="not pathogen"
			  -> "S"    <LIST>: A list with 7 items:
					[mean <FLOAT>, sd <FLOAT>, mc_error <FLOAT>, hpd_2.5 <FLOAT>, hpd_97.5 <FLOAT>, p_lt_0 <FLOAT>, p_gt_0 <FLOAT>]
			  -> "SL"   <INT> : genome sequence length
			  -> [rank] <DICT>: # LCR rank
		   			-> "R"     <INT>  : raw read count
					-> "M"     <INT>  : mapped length (match + mismatch)
					-> "Mnr"   <INT>  : mapped length / num_refs
					-> "T"     <INT>  : total read length
					-> "N"     <INT>  : total mismatch bases
					-> "Rnb"   <FLOAT>: normalize read count by both identity and # of ref (1*identity/num_refs)
					-> "Rnr"   <FLOAT>: normalize read count by the number of references (1/num_refs)
					-> "Ri"    <FLOAT>: identity (mapped_length-nm)/read_length

	  res_tree <DICT>:
		  -> pid <DICT> : parent taxid
			  -> tid <STR> : tid

	  res_rank <DICT>:
		  -> rank <DICT>: rank
			  -> "t_rsMR"   <INT>   : total rank-specific read count
			  -> "t_rsMRNr" <FLOAT> : total rank-specific RNR
			  -> "t_nsMR"   <INT>   : total non-rank-specific read count
			  -> "t_nsMRNr" <FLOAT> : total non-rank-specific RNR
	"""
	res_rollup = t._autoVivification() 
	res_rollup_str = t._autoVivification()
	res_rank = t._autoVivification()
	res_tree = t._autoVivification()
	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom','root']
	str_rollup_cnt = 0

	# rollup to strain first
	for ref in r:
		(acc, length, tid, flag) = ref.split('|')
		rank = t.taxid2rank(tid)
		if rank != "strain": continue

		if tid in res_rollup_str:
			#res_rollup_str[tid]["ML"] += ";%s:%s" % (ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
			res_rollup_str[tid]["MB"] += r[ref]["MB"]
			res_rollup_str[tid]["TB"] += r[ref]["TB"]
			res_rollup_str[tid]["MR"] += r[ref]["MR"]
			res_rollup_str[tid]["PR"] += r[ref]["PR"]
			res_rollup_str[tid]["NM"] += r[ref]["NM"]
			res_rollup_str[tid]["LL"] += r[ref]["LL"]
		else:
			#res_rollup_str[tid]["ML"] = "%s:%s" % (ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
			res_rollup_str[tid]["MB"] = r[ref]["MB"]
			res_rollup_str[tid]["TB"] = r[ref]["TB"]
			res_rollup_str[tid]["MR"] = r[ref]["MR"]
			res_rollup_str[tid]["PR"] = r[ref]["PR"]
			res_rollup_str[tid]["NM"] = r[ref]["NM"]
			res_rollup_str[tid]["LL"] = r[ref]["LL"]
			res_rollup_str[tid]["F"]  = flag[0]

		for lcr_lvl in ranks:
			if lcr_lvl in r[ref]:
				if not lcr_lvl in res_rollup_str[tid]:
					res_rollup_str[tid][lcr_lvl]["R"]     = 0 
					res_rollup_str[tid][lcr_lvl]["M"]     = 0 
					res_rollup_str[tid][lcr_lvl]["Mnr"]   = 0 
					res_rollup_str[tid][lcr_lvl]["T"]     = 0 
					res_rollup_str[tid][lcr_lvl]["N"]     = 0 
					res_rollup_str[tid][lcr_lvl]["Rnb"]   = 0 
					res_rollup_str[tid][lcr_lvl]["Rnr"]   = 0 
				
				res_rollup_str[tid][lcr_lvl]["R"]     += r[ref][lcr_lvl]["R"]
				res_rollup_str[tid][lcr_lvl]["M"]     += r[ref][lcr_lvl]["M"]
				res_rollup_str[tid][lcr_lvl]["Mnr"]   += r[ref][lcr_lvl]["Mnr"]
				res_rollup_str[tid][lcr_lvl]["T"]     += r[ref][lcr_lvl]["T"]
				res_rollup_str[tid][lcr_lvl]["N"]     += r[ref][lcr_lvl]["N"]
				res_rollup_str[tid][lcr_lvl]["Rnb"]   += r[ref][lcr_lvl]["Rnb"]
				res_rollup_str[tid][lcr_lvl]["Rnr"]   += r[ref][lcr_lvl]["Rnr"]
	
	print_message( "%s strain(s) mapped." % len(res_rollup_str), argvs.silent, begin_t, logfile )

	# apply cutoffs strain level and rollup to higher levels
	for taxid in res_rollup_str:
		if not taxid in genome_size:
			genome_size[taxid]

		# calculate strain DC (depth of coverage)
		res_rollup_str[taxid]["DC"] = res_rollup_str[taxid]["MB"]/genome_size[taxid]
		res_rollup_str[taxid]["DCnr"] = 0
		
		# calculate strain DCnr and DC at each rank
		for lcr_lvl in ranks:
			if lcr_lvl in res_rollup_str[taxid]:
				res_rollup_str[taxid][lcr_lvl]["DC"] = res_rollup_str[taxid][lcr_lvl]["M"]/genome_size[taxid]
				res_rollup_str[taxid][lcr_lvl]["DCnr"] = res_rollup_str[taxid][lcr_lvl]["Mnr"]/genome_size[taxid]
				res_rollup_str[taxid]["DCnr"] += res_rollup_str[taxid][lcr_lvl]["DCnr"]

		# calculate strain linear coverage
		res_rollup_str[taxid]["LC"] = res_rollup_str[taxid]["LL"]/genome_size[taxid]

		# pass genome size
		res_rollup_str[taxid]["SL"] = genome_size[taxid]

		# filtering out non-qualified strains
		#if mc > res_rollup_str[taxid]["LC"] or \
		#   mr > res_rollup_str[taxid]["MR"] or \
		#   ml > res_rollup_str[taxid]["LL"] or \
		#   md > res_rollup_str[taxid]["DC"]:
		#   #(mb > res_rollup_str[taxid]["strain"]["Rnb"] if "strain" in res_rollup_str[taxid] else -1 ) or \
		#   continue

		#tree = t.taxid2fullLinkDict( taxid )
		lng = t.taxid2lineageDICT(taxid, 1, 1)
		str_rollup_cnt += 1

		# rollup to each ranks (including strain)
		for rank in lng:
			if rank == "type": continue #skip "type" level
			#res_tree[pid][tid] = 1
			tid = lng[rank]['taxid']

			# if no particular rank for this organism, use strain
			if tid == 0:
				tid = lng[rank]['name'] # for example: 'Flaviviridae - no_o_rank - no_c_rank'

			if tid in res_rollup:
				#res_rollup[tid]["ML"] += ";%s" % res_rollup_str[taxid]["ML"]
				res_rollup[tid]["MB"] += res_rollup_str[taxid]["MB"]
				res_rollup[tid]["TB"] += res_rollup_str[taxid]["TB"]
				res_rollup[tid]["MR"] += res_rollup_str[taxid]["MR"]
				res_rollup[tid]["PR"] += res_rollup_str[taxid]["PR"]
				res_rollup[tid]["NM"] += res_rollup_str[taxid]["NM"]
				res_rollup[tid]["LL"] += res_rollup_str[taxid]["LL"]
				res_rollup[tid]["DC"] += res_rollup_str[taxid]["DC"]
				res_rollup[tid]["DCnr"] += res_rollup_str[taxid]["DCnr"]
				res_rollup[tid]["SL"] += res_rollup_str[taxid]["SL"]
				if res_rollup_str[taxid]["LC"] > res_rollup[tid]["LC"]:
					res_rollup[tid]["LC"] = res_rollup_str[taxid]["LC"] 
			else:
				#res_rollup[tid]["ML"] = res_rollup_str[taxid]["ML"]
				res_rollup[tid]["LVL"] = rank
				res_rollup[tid]["MB"] = res_rollup_str[taxid]["MB"]
				res_rollup[tid]["TB"] = res_rollup_str[taxid]["TB"]
				res_rollup[tid]["MR"] = res_rollup_str[taxid]["MR"]
				res_rollup[tid]["PR"] = res_rollup_str[taxid]["PR"]
				res_rollup[tid]["NM"] = res_rollup_str[taxid]["NM"]
				res_rollup[tid]["LL"] = res_rollup_str[taxid]["LL"]
				res_rollup[tid]["DC"] = res_rollup_str[taxid]["DC"]
				res_rollup[tid]["DCnr"] = res_rollup_str[taxid]["DCnr"]
				res_rollup[tid]["SL"] = res_rollup_str[taxid]["SL"]
				res_rollup[tid]["LC"] = res_rollup_str[taxid]["LC"]
				res_rollup[tid]["F"]  = res_rollup_str[taxid]["F"]

			for lcr_lvl in ranks:
				if lcr_lvl in res_rollup_str[taxid]:
					if not lcr_lvl in res_rollup[tid]:
						res_rollup[tid][lcr_lvl]["R"]     = 0 
						res_rollup[tid][lcr_lvl]["M"]     = 0 
						res_rollup[tid][lcr_lvl]["Mnr"]   = 0 
						res_rollup[tid][lcr_lvl]["T"]     = 0 
						res_rollup[tid][lcr_lvl]["N"]     = 0
						res_rollup[tid][lcr_lvl]["DC"]    = 0 
						res_rollup[tid][lcr_lvl]["DCnr"]  = 0 
						res_rollup[tid][lcr_lvl]["Rnb"]   = 0 
						res_rollup[tid][lcr_lvl]["Rnr"]   = 0 
					
					res_rollup[tid][lcr_lvl]["R"]     += res_rollup_str[taxid][lcr_lvl]["R"]
					res_rollup[tid][lcr_lvl]["M"]     += res_rollup_str[taxid][lcr_lvl]["M"]
					res_rollup[tid][lcr_lvl]["Mnr"]   += res_rollup_str[taxid][lcr_lvl]["Mnr"]
					res_rollup[tid][lcr_lvl]["T"]     += res_rollup_str[taxid][lcr_lvl]["T"]
					res_rollup[tid][lcr_lvl]["N"]     += res_rollup_str[taxid][lcr_lvl]["N"]
					res_rollup[tid][lcr_lvl]["DC"]    += res_rollup_str[taxid][lcr_lvl]["DC"]
					res_rollup[tid][lcr_lvl]["DCnr"]  += res_rollup_str[taxid][lcr_lvl]["DCnr"]
					res_rollup[tid][lcr_lvl]["Rnb"]   += res_rollup_str[taxid][lcr_lvl]["Rnb"]
					res_rollup[tid][lcr_lvl]["Rnr"]   += res_rollup_str[taxid][lcr_lvl]["Rnr"]

	#print_message( "%s strain(s) rolled up."%str_rollup_cnt, argvs.silent, begin_t, logfile )

	# Generate following values:
	# 1) normalized read count
	# 2) RPKM
	# 3) check and add pathogen info

	for tid in res_rollup:
		trank =  res_rollup[tid]["LVL"]
		lvl = ranks.index(trank) if trank in ranks else -1

		# init
		res_rollup[tid]["rsMRNb"] = 0
		res_rollup[tid]["MRNr"]   = 0
		res_rollup[tid]["rsMR"]   = 0
		res_rollup[tid]["rsMB"]   = 0
		res_rollup[tid]["rsMRNr"] = 0
		res_rollup[tid]['rsDCnr'] = 0
		res_rollup[tid]['rsDC'] = 0
		res_rollup[tid]["nsMR"]   = 0
		res_rollup[tid]["nsMRNr"] = 0

		if trank in ranks:
			if not trank in res_rank:
				res_rank[trank]["t_rsMR"]   = 0
				res_rank[trank]["t_rsMRNr"] = 0 
				res_rank[trank]["t_nsMR"]   = 0
				res_rank[trank]["t_nsMRNr"] = 0

		# add up read counts and mapped bases normalized by # of references
		for lcr_lvl in ranks:
			if lcr_lvl in res_rollup[tid]:
				# add up read counts and mapped bases normalized by # of references
				res_rollup[tid]["MRNr"] += res_rollup[tid][lcr_lvl]["Rnr"]
				# calculate read-mapping identity
				res_rollup[tid][lcr_lvl]["Ri"] = (res_rollup[tid][lcr_lvl]["M"]-res_rollup[tid][lcr_lvl]["N"])/res_rollup[tid][lcr_lvl]["T"]

		# only add up rank specific normalized read count for major ranks
		rs_flag = 1
		if trank in ranks:
			for idx, lcr_lvl in enumerate(ranks):
				# set rank specific flag
				if idx > lvl: rs_flag = 0

				if lcr_lvl in res_rollup[tid]:
					if rs_flag:
						res_rollup[tid]['rsMRNb']   += res_rollup[tid][lcr_lvl]["Rnb"]
						res_rollup[tid]['rsMR']     += res_rollup[tid][lcr_lvl]["R"]
						res_rollup[tid]['rsMB']     += res_rollup[tid][lcr_lvl]["M"]
						res_rollup[tid]['rsMRNr']   += res_rollup[tid][lcr_lvl]["Rnr"]
						res_rollup[tid]['rsDCnr']   += res_rollup[tid][lcr_lvl]["DCnr"]
						res_rollup[tid]['rsDC']     += res_rollup[tid][lcr_lvl]["DC"]
						res_rank[trank]['t_rsMR']   += res_rollup[tid][lcr_lvl]["R"]
						res_rank[trank]['t_rsMRNr'] += res_rollup[tid][lcr_lvl]["Rnr"]
					else:
						res_rollup[tid]['nsMR']     += res_rollup[tid][lcr_lvl]["R"]
						res_rollup[tid]['nsMRNr']   += res_rollup[tid][lcr_lvl]["Rnr"]
						res_rank[trank]['t_nsMR']   += res_rollup[tid][lcr_lvl]["R"]
						res_rank[trank]['t_nsMRNr'] += res_rollup[tid][lcr_lvl]["Rnr"]

		# calculating RPKM
		# In order to keep the sum of read count equivalent to total number of real mapped reads, 
		# MRNr is used to calculate RPKM
		try:
			res_rollup[tid]["RPKM"] = res_rollup[tid]["MRNr"] / (res_rollup[tid]["SL"]/1e3) / (tol_num_mapped_reads/1e6)
		except:
			res_rollup[tid]["RPKM"] = 0

		# check pathogen and add info
		if patho_meta:
			if tid in patho_meta:
				res_rollup[tid]["P"] = tid
			elif tid.split('.')[0] in patho_meta:
				res_rollup[tid]["P"] = tid.split('.')[0]
			elif trank == "strain" and t.taxid2taxidOnRank(tid, "species") in patho_meta:
				res_rollup[tid]["P"] = t.taxid2taxidOnRank(tid, "species")
			else:
				res_rollup[tid]["P"] = 0
	return res_rollup

def scoreStandalone(res_rollup, uniq_meta, method):
	"""scoreStandalone
	Calculate score using rank uniqueness data and add the score to res_rollup
	"""

	for cur_tid in res_rollup:
		trank = res_rollup[cur_tid]["LVL"]
		if not "S_SA"    in res_rollup[cur_tid]: res_rollup[cur_tid]["S_SA"] = 0

		cur_score = s.scoreStrainUniqueness(res_rollup[cur_tid], cur_tid, uniq_meta, method)
		res_rollup[cur_tid]["S_SA_CL"] = cur_score
	
		if not "-" in cur_tid:
			lng = t.taxid2lineageDICT(cur_tid, 1, 1)

			# if the score of the 
			for rank in lng:
				if rank=="type": continue
				#res_tree[pid][tid] = 1
				tid = lng[rank]['taxid']
				
				if tid == 0:
					tid = lng[rank]['name'] # for example: 'Flaviviridae - no_o_rank - no_c_rank'

				if not tid in res_rollup: continue
				
				if not "S_SA" in res_rollup[tid]:
					res_rollup[tid]["S_SA"] = cur_score
				else:
					if res_rollup[tid]["S_SA"] == "NA":
						res_rollup[tid]["S_SA"] = cur_score
					elif cur_score == "NA":
						continue
					# if the genome coverage contributed strain has a higer score, replace the current score
					elif res_rollup[tid]["S_SA"] < cur_score and res_rollup[tid]["LC"] <= res_rollup[cur_tid]["LC"]:
						res_rollup[tid]["S_SA"] = cur_score

	return res_rollup

def outputResultsAsReport(res_rollup, o, relAbu, user_taxid, display_fields, score_method, ms, mc, mb, mr, ml, md, mrd, display_all):
	"""outputResultsAsReport

	Arguments:
	  Arg1:
	  res_rollup <DICT>: read mapping results roll up to superkingdom
		  -> taxid <DICT>: taxid as keys
			  -> "ML"   <LIST>: list mapped reference(s) and their locations
			  -> "MB"   <INT> : # of mapped bases (mapping region only -> match + mismatch)
			  -> "TB"   <INT> : total length of mapped reads (match + mismatch)
			  -> "MR"   <INT> : # of mapped reads
			  -> "PR"   <INT> : # of primary mapped reads
			  -> "DC"   <FLOAT> : depth of coverage
			  -> "DCnr" <FLOAT> : depth of coverage normalized by the number of references
			  -> "rsDCnr" <FLOAT> : rank specific depth of coverage normalized by the number of references
			  -> "MRNr" <INT> : Total normalize read count by the number of references
			  -> "rsMRNb" <INT> : Total rank specific normalized reads
			  -> "rsMB" <INT> : Total rank specific mapped bases
			  -> "RPKM" <INT> : RPKM
			  -> "NM"   <INT> : # of mismatches
			  -> "LL"   <INT> : linear length (bp)
			  -> "LC"   <INT> : best linear coverage (prop)
			  -> "P"    <STR> : Pathogen taxonomy id; 0="not pathogen"
			  -> "S"    <LIST>: A list with 7 items:
					[mean <FLOAT>, sd <FLOAT>, mc_error <FLOAT>, hpd_2.5 <FLOAT>, hpd_97.5 <FLOAT>, p_lt_0 <FLOAT>, p_gt_0 <FLOAT>]
			  -> "SL"   <INT> : genome sequence length
			  -> [rank] <DICT>: # LCR rank
		   			-> "R"     <INT>  : raw read count
					-> "M"     <INT>  : mapped length (match + mismatch)
					-> "Mnr"   <INT>  : mapped length / num_refs
					-> "T"     <INT>  : total read length
					-> "N"     <INT>  : total mismatch bases
					-> "Rnb"   <FLOAT>: normalize read count by both identity and # of ref (Rni/num_refs)
					-> "Rnr"   <FLOAT>: normalize read count by the number of references (1/num_refs)
					-> "Rni"   <FLOAT>: normalize read count by identity (mapped_length-nm)/read_length
					-> "R_IDT" <LIST> : list of the identity of each read
					-> "R_MAT" <LIST> : list of the # of matches of each read
					-> "R_SCR" <LIST> : list of alignment scores
	  Arg2: o <OBJ>
	  Arg3: relAbu <DICT>
	  Arg4: display_fields <LIST>
	  Arg5: ms <FLOAT>
	  Arg6: mr <FLOAT>
	  Arg7: ml <FLOAT>
	
	Return: None
	"""
	output = t._autoVivification()
	note = ""
	major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}
	lcr_ranks = ['strain','species','genus','family','order','class','phylum','superkingdom','root']

	#convert res_rollup to an output dict that maps keys to the rank
	for tid in res_rollup:
		#skip 'cellular organisms', it's mistakly reported as 'strain'
		if tid == "131567": continue

		# user_taxid
		rank = res_rollup[tid]["LVL"]

		# get score
		score = "NA"
		score_bg = "none"
		score_sa = "none"
		if 'S_BG' in res_rollup[tid]:
			score_bg = res_rollup[tid]["S_BG"]
		if 'S_SA' in res_rollup[tid]:
			score_sa = res_rollup[tid]["S_SA"]

		if score_method == "bg" and score_bg != "none":
			score = score_bg
		elif score_method == "standalone" and score_sa != "none":
			score = score_sa
		else:
			if type(score_bg) is str: score_bg = 0
			if type(score_sa) is str: score_sa = 0
			score = score_sa * score_bg

		# calculate abundance for each taxonomy and total abundance for each rank 
		if rank in major_ranks:
			abundance=0
			if (str(score) != "NA" and ms > float(score) or
					mb > res_rollup[tid]["rsMRNb"] or
					mr > res_rollup[tid]["MR"] or
					ml > res_rollup[tid]["LL"] or
					md > res_rollup[tid]["DC"] or
					mc > res_rollup[tid]["LC"] or
					mrd > res_rollup[tid]["rsDCnr"] ):
				abundance = 0
			else:
				if relAbu == "TOTAL_BP_MAPPED":
					abundance = res_rollup[tid]["MB"]
				elif relAbu == "READ_COUNT":
					abundance = res_rollup[tid]["MR"]
				elif relAbu == "READ_COUNT_RNR":
					abundance = res_rollup[tid]["MRNr"]
				elif relAbu == "READ_COUNT_RSNB":
					abundance = res_rollup[tid]["rsMRNb"]
				elif relAbu == "PRI_READ_COUNT":
					abundance = res_rollup[tid]["PR"]
				elif relAbu == "RPKM":
					abundance = res_rollup[tid]["RPKM"]
				else:
					abundance = res_rollup[tid]["rsDCnr"]

			output[rank]["RES"][tid] = abundance

			if "TOT_ABU" in output[rank]:
				output[rank]["TOT_ABU"] += abundance
			else:
				output[rank]["TOT_ABU"] = abundance

	# prepare report headers
	# report filed sets: ['basic', 'r', rnb', 'rnr', 'ri', 'patho', 'score', 'ref', 'full', 'all']
	basic_fields = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"LEVEL", "NAME", "TAXID", "READ_COUNT", "READ_COUNT_RNR", 
						"READ_COUNT_RSNB", "LINEAR_COV", "DEPTH_COV", "DEPTH_COV_NR", "RS_DEPTH_COV_NR",
						"PATHOGEN", "SCORE", "REL_ABUNDANCE" )

	full_fields  = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"ABUNDANCE", "TOTAL_BP_MISMATCH", "NOTE", "RPKM", "PRI_READ_COUNT",
						"TOL_RS_READ_CNT", "TOL_NS_READ_CNT", "TOL_RS_RNR", "TOL_NS_RNR", "TOL_GENOME_SIZE",
						"LINEAR_LENGTH", "TOTAL_BP_MAPPED", "RS_DEPTH_COV", "FLAG"
						) if "full" in display_fields else ""
	
	# rank specific data
	r_fields   = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"STR", "SPE", "GEN", "FAM", "ORD", "CLA", "PHY", "SK", "ROOT" ) if "r" in display_fields else ""
	rnb_fields = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"STR_rnb", "SPE_rnb", "GEN_rnb", "FAM_rnb", "ORD_rnb", "CLA_rnb", "PHY_rnb", "SK_rnb", "ROOT_rnb" ) if "rnb" in display_fields else ""
	rnr_fields = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"STR_rnr", "SPE_rnr", "GEN_rnr", "FAM_rnr", "ORD_rnr", "CLA_rnr", "PHY_rnr", "SK_rnr", "ROOT_rnr" ) if "rnr" in display_fields else ""
	ri_fields  = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						"STR_ri", "SPE_ri", "GEN_ri", "FAM_ri", "ORD_ri", "CLA_ri", "PHY_ri", "SK_ri", "ROOT_ri" ) if "ri" in display_fields else ""

	patho_fields  = "\tSOURCE\tLOCATION\tHOST\tDISEASE" if "patho" in display_fields else ""

	score_fields  = "\tSCORE_UNIQ\tSCORE_BG\tSCORE_UNIQ_CUR_LVL" if "score" in display_fields else ""

	# output report headers
	o.write( "%s%s%s%s%s%s%s%s\n" % (
		basic_fields, full_fields, r_fields, rnb_fields, rnr_fields, ri_fields, patho_fields, score_fields))

	# prepare field values
	for rank in sorted( major_ranks, key=major_ranks.__getitem__ ):
		
		taxas = output[rank]["RES"] #all corresponding abundance in this rank

		for tid in sorted( taxas, key=taxas.__getitem__, reverse=True ):

			# check if any ref is belong to provided taxid
			if user_taxid and not isDescendant(tid, user_taxid): continue

			# get score
			score = "NA"
			score_bg = "none"
			score_sa = "none"
			if 'S_BG' in res_rollup[tid]:
				score_bg = res_rollup[tid]["S_BG"]
			if 'S_SA' in res_rollup[tid]:
				score_sa = res_rollup[tid]["S_SA"]

			if score_method == "bg" and score_bg != "none":
				score = score_bg
			elif score_method == "standalone" and score_sa != "none":
				score = score_sa
			else:
				if type(score_bg) is str: score_bg = 0
				if type(score_sa) is str: score_sa = 0
				score = score_sa * score_bg

			note = ""
			note += "filtered out by minScore. "  if str(score) != "NA" and ms > float(score) else ""
			note += "filtered out by minCov. "    if mc > res_rollup[tid]["LC"] else ""
			note += "filtered out by minRsnb. "   if mb > res_rollup[tid]["rsMRNb"] else ""
			note += "filtered out by minReads. "  if mr > res_rollup[tid]["MR"] else ""
			note += "filtered out by minLen. "    if ml > res_rollup[tid]["LL"] else ""
			note += "filtered out by minDc. "     if md > res_rollup[tid]["DC"] else ""
			note += "filtered out by minRsdcnr. " if mrd > res_rollup[tid]["rsDCnr"] else ""

			if note and not display_all: continue

			#basic_fields, full_fields, r_fields, rnb_fields, rnr_fields, ri_fields, patho_fields, score_fields, ref_fields

			# BASIX FIELDS
			# values for basic fields
			basic_fields = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
				rank,                                                    # 1. "LEVEL"        
				t.taxid2name(tid) if not " - no_" in tid else tid,       # 2. "NAME"                     
				tid if not " - no_" in tid else 0,                       # 3. "TAXID"                                        
				res_rollup[tid]["MR"],                                   # 4. "READ_COUNT"                         
				"%.2f"%res_rollup[tid]["MRNr"],                          # 5. "READ_COUNT_RNR"                                  
				"%.2f"%res_rollup[tid]["rsMRNb"],                        # 6. "READ_COUNT_RSNB"                                    
				"%.4f"%res_rollup[tid]["LC"],                            # 7. "LINEAR_COV"      
				"%.4f"%res_rollup[tid]["DC"],                            # 8. "DEPTH_COV"                          
				"%.4f"%res_rollup[tid]["DCnr"],                          # 9. "DEPTH_COV_NR"                                  
				"%.4f"%res_rollup[tid]["rsDCnr"],                        # 10."RS_DEPTH_COV_NR"                                
				"Pathogen" if res_rollup[tid]["P"] else "",              # 11."PATHOGEN"                                              
				"%.4f"%score if score != "NA" else "NA",                 # 12."SCORE"                                           
				0 if (note or not output[rank]["TOT_ABU"]) else "%.4f"%(output[rank]["RES"][tid]/output[rank]["TOT_ABU"]) # 13."REL_ABUNDANCE"                                                            
			)

			# FULL FIELDS
			# values for full_fields
			full_fields = "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
				output[rank]["RES"][tid],       # 14."ABUNDANCE"               
				res_rollup[tid]["NM"],          # 15."TOTAL_BP_MISMATCH"            
				note,                           # 16."NOTE"    
				"%.4f"%res_rollup[tid]["RPKM"], # 17."RPKM"                     
				res_rollup[tid]["PR"],          # 18."PRI_READ_COUNT"                         
				res_rollup[tid]["rsMR"],        # 19."TOL_RS_READ_CNT"              
				res_rollup[tid]["nsMR"],        # 20."TOL_NS_READ_CNT"              
				res_rollup[tid]["rsMRNr"],      # 21."TOL_RS_RNR",                
				res_rollup[tid]["nsMRNr"],      # 22."TOL_NS_RNR"
				res_rollup[tid]["SL"],          # 23."TOL_GENOME_SIZE",                
				res_rollup[tid]["LL"],          # 24."LINEAR_LENGTH"                
				res_rollup[tid]["MB"],          # 25."TOTAL_BP_MAPPED"              
				res_rollup[tid]["rsDC"],        # 26."RS_DEPTH_COV"    
				res_rollup[tid]["F"],           # 27."FLAG"              
			) if "full" in display_fields else ""

			#prepare lcr fields
			#r_fields, rnb_fields, rnr_fields, ri_fields
			r_fields   = ""
			rnb_fields = ""
			rnr_fields = ""
			ri_fields = ""
			for lcr_lvl in lcr_ranks:
				if lcr_lvl in res_rollup[tid]:
					r_fields   += "\t%s"   % ( res_rollup[tid][lcr_lvl]["R"]   ) if "r" in display_fields else ""
					rnb_fields += "\t%.2f" % ( res_rollup[tid][lcr_lvl]["Rnb"] ) if "rnb" in display_fields else ""
					rnr_fields += "\t%.2f" % ( res_rollup[tid][lcr_lvl]["Rnr"] ) if "rnr" in display_fields else ""
					ri_fields  += "\t%.4f" % ( res_rollup[tid][lcr_lvl]["Ri"] ) if "ri" in display_fields else ""
				else:
					r_fields   += "\t0" if "r" in display_fields else ""
					rnb_fields += "\t0" if "rnb" in display_fields else ""
					rnr_fields += "\t0" if "rnr" in display_fields else ""
					ri_fields  += "\t0" if "ri" in display_fields else ""

			# pathogen
			patho_fields = ""
			if res_rollup[tid]["P"]:
				p_tid = res_rollup[tid]["P"]
				patho_fields  = "\t%s\t%s\t%s\t%s" % (
					patho_meta[p_tid]['source'],
					patho_meta[p_tid]['location'],
					patho_meta[p_tid]['host'],
					patho_meta[p_tid]['disease']
				) if "patho" in display_fields else ""
			else:
				patho_fields  = "\t\t\t\t" if "patho" in display_fields else ""

			# score
			score_fields = ""
			if "S_SA" in res_rollup[tid]:
				s = "%.6f"%res_rollup[tid]["S_SA"] if res_rollup[tid]["S_SA"] != "none" else "none"
				score_fields += "\t%s"%s if "score" in display_fields else ""
			else:
				score_fields += "\t" if "score" in display_fields else ""

			if "S_BG" in res_rollup[tid]:
				s = "%.6f"%res_rollup[tid]["S_BG"] if res_rollup[tid]["S_BG"] != "none" else "none"
				score_fields += "\t%s"%s if "score" in display_fields else ""
			else:
				score_fields += "\t" if "score" in display_fields else ""

			if "S_SA_CL" in res_rollup[tid]:
				s = "%.6f"%res_rollup[tid]["S_SA_CL"] if res_rollup[tid]["S_SA_CL"] != "none" else "none"
				score_fields += "\t%s"%s if "score" in display_fields else ""
			else:
				score_fields += "\t" if "score" in display_fields else ""

			# output result
			o.write( "%s%s%s%s%s%s%s%s\n" % (
				basic_fields, full_fields, r_fields, rnb_fields, rnr_fields, ri_fields, patho_fields, score_fields))

def outputResultsAsLineage(res_rollup, o, relAbu, mode, score_method, ms, mb, mr, ml, md, display_all):
	for tid in res_rollup:
		rank = res_rollup[tid]["LVL"]
		if rank != "species": continue
		if " - no_" in tid: continue

		# get score
		score = "NA"
		score_bg = "NA"
		score_sa = "NA"

		if "S_BG" in res_rollup[tid] and type(res_rollup[tid]["S_BG"]) is float:
			score_bg = res_rollup[tid]["S_BG"]

		if "S_SA" in res_rollup[tid] and type(res_rollup[tid]["S_SA"]) is float:
			score_sa = res_rollup[tid]["S_SA"]

		if score_method == "bg":
			score = score_bg
		elif score_method == "standalone":
			score = score_sa
		else:
			if score_bg != "NA" and score_sa != "NA":
				score = score_sa * score_bg

		if not display_all and not score == "NA" and ms > score: continue
		if not display_all and mb > res_rollup[tid]["rsMRNb"]: continue
		if not display_all and mr > res_rollup[tid]["MR"]: continue
		if not display_all and ml > res_rollup[tid]["LL"]: continue
		if not display_all and md > res_rollup[tid]["DC"]: continue

		if relAbu == "LINEAR_LENGTH":
			abundance = res_rollup[tid]["LL"]
		elif relAbu == "TOTAL_BP_MAPPED":
			abundance = res_rollup[tid]["MB"]
		elif relAbu == "READ_COUNT":
			abundance = res_rollup[tid]["MR"]
		elif relAbu == "PRI_READ_COUNT":
			abundance = res_rollup[tid]["PR"]
		elif relAbu == "READ_COUNT_RNR":
			abundance = res_rollup[tid]["MRNr"]
		elif relAbu == "READ_COUNT_RSNB":
			abundance = res_rollup[tid]["rsMRNb"]
		elif relAbu == "RPKM":
			abundance = res_rollup[tid]["RPKM"]
		else:
			abundance = res_rollup[tid]["DC"]

		lineage = t.taxid2lineage(tid)
		o.write( "%s\t%s\n" %
			( abundance,
				'\t'.join( lineage.split('|') )
		))

def mergingSAM( in_sam_files, out_sam_file, min_score=0, host_tag="H", nm_range=0 ):
	"""
	This function is used to merge different SAM files.

	ARGUMENT:
		inSamFin_sam_files LIST   filename(s) of input SAM files
		outSout_sam_file   OBJ    file object for output merged SAM files
		min_score          INT
		host_tag           STR
		ignore_tag         LIST
		nm_range           INT
	RETURN:
		tol_host           INT    total host reads
		tol_mapped         INT    total mapped reads
	"""
	sam={}
	tol_host=0
	tol_rootspec=0
	tol_mapped=0

	def purifyXA( xa, nm, nm_range):
		"""
		Clean up XAs that do not pass nm_range
		"""
		xa_str = xa.replace("XA:Z:","")
		idv_xas = xa_str.rstrip(";").split(";")
		tmp = list(idv_xas)

		for idv_xa in tmp:
			if not idv_xa.endswith( ","+str(nm) ):
				idv_xa_nm = idv_xa.split(",")[-1]
				if int(idv_xa_nm) > nm+nm_range:
					idv_xas.remove( idv_xa )

		if len(idv_xas):
			return "XA:Z:"+";".join(idv_xas)+";"
		else:
			return ""

	def prepFinalSAM( sam, readname, mate, nm_range):
		"""
		Prepare final SAM record
		"""
		aln = sam[readname][mate]["SAM"]
		# purify XA tag
		if "XA:Z" in aln:
			tmp = aln.split("\t")
			new_xa = purifyXA( tmp[-1], sam[readname][mate]["NM"], nm_range )
			# if XA tag exists after purification
			if new_xa:
				tmp[-1]=new_xa
			else:
				tmp = tmp[:-1]
			# update record
			sam[readname][mate]["SAM"] = "\t".join(tmp)
			aln = sam[readname][mate]["SAM"]
		
		# filter out HOST reads
		if "|"+host_tag in aln:
			return "HOST"

		# filter out ROOT-specific reads
		if "XA:Z" in aln:
			# a read can map to more than 1 superkingdom (category) is a ROOT-specific read
			all_sk = set(findall('\|(\w)\w*[,\t]', aln))
			if len(all_sk)>1:
				return "ROOTSPEC"

		# output sam record
		return sam[readname][mate]["SAM"]

	for filename in in_sam_files:
		with open( filename ) as f:
			for line in f:
				"""
				Example input:
				@SQ     SN:ref1 LN:40
				@SQ     SN:ref2 LN:40
				read1	0	ref2	11	0	10M	*	0	0	GGGGGGGGGG	*	NM:i:0	MD:Z:10	AS:i:10	XS:i:10	XA:Z:ref3,+21,10M,0;ref1,+1,10M,0;ref4,+31,10M,0;
				read2	16	ref4	31	0	10M	*	0	0	GGGGGGGGGG	*	NM:i:0	MD:Z:10	AS:i:10	XS:i:10	XA:Z:ref1,-1,10M,0;ref3,-21,10M,0;ref2,-11,10M,0;	
				XA tag contains suboptional alignments in (chr,pos,CIGAR,NM;) format
				"""

				line = line.strip()
				temp = line.split('\t')
				samflag = int(temp[1])
				score = -1
				readname = temp[0]
				strand = "-" if samflag & 16 else "+"

				# auto detect paired-end or single-end reads
				if argvs.singleEnd == "auto":
					if int(temp[1]) & 1:
						argvs.singleEnd = False
					else:
						argvs.singleEnd = True

				# skip paired-end reads but flagged mate unmapped (0x8)
				if not argvs.singleEnd and samflag & 8:
					continue

				# add seq prefix to the readname if the read is single-end
				if not samflag & 1:
					seq = temp[9]
					if strand == "-":
						seq = seqReverseComplement(temp[9])
					# paired-end reads but forced to treat as signle-end reads
					if int(temp[1])&1 and argvs.singleEnd:
						readname += seq[0:5] # add the first 5bp of its sequence to separate
					
				# find out mate number
				mate = "1"
				if samflag & 128:
					mate = "2"

				# get alignment score
				if "AS:i" in line:
					s_as = search('AS:i:(\d+)', line)
					s_ms = search('ms:i:(\d+)', line)
					try:
						score = int(s_as.group(1))
					except:
						try:
							score = int(s_ms.group(1))
						except:
							score = 0

				# skip if alignment score doesn't reach cutoff
				if score < min_score:
					continue

				# get NM
				s_nm = search('NM:i:(\d+)', line)
				nm   = int(s_nm.group(1))
		
				if readname in sam and mate in sam[readname]:
					# when new best score found, replace current record
					if score > sam[readname][mate]["SCR"]:
						sam[readname][mate]["SCR"] = score
						sam[readname][mate]["SAM"] = line
						sam[readname][mate]["SRD"] = strand
						sam[readname][mate]["NM"]  = nm
					# equal best score found, add current alignments to existing record
					elif score == sam[readname][mate]["SCR"]:
						# get strand
						srd = "-" if samflag & 16 else "+"
						# get XA
						xa = ""
						if "XA:Z" in line:
							s_xa = search('XA:Z:(\S+)', line)
							xa = s_xa.group(1)
						# convert primary alignment to XA
						pa = "%s,%s%s,%s,%s;"%(temp[2],srd,temp[3],temp[5],nm)
						xa = xa + pa
						# change strand if the primary alignment different from the record
						if srd != sam[readname][mate]["SRD"]:
							xa = xa.replace(",+", ",^")
							xa = xa.replace(",-", ",+")
							xa = xa.replace(",^", ",-")
						# merge to existing record
						if "XA:Z" in sam[readname][mate]["SAM"]:
							sam[readname][mate]["SAM"] += xa
						else:
							sam[readname][mate]["SAM"] += "\tXA:Z:%s"%xa
					else:
						pass
				else:
					if not readname in sam:
						sam[readname] = {}
					if not mate in sam[readname]:
						sam[readname][mate] = {}

					sam[readname][mate]["SCR"] = score
					sam[readname][mate]["SAM"] = line
					sam[readname][mate]["SRD"] = strand
					sam[readname][mate]["NM"]  = nm
	
	# save merged alignments
	for readname in sam:
		for mate in sam[readname]:
			tol_mapped += 1
			final_sam = prepFinalSAM(sam, readname, mate, nm_range)
			
			if final_sam == "HOST":
				# filter out HOST reads
				tol_host += 1			
			elif final_sam == "ROOTSPEC":
				# filter out ROOT-specific reads
				tol_rootspec += 1
			else:
				# output sam record
				out_sam_file.write("%s\n"%final_sam)

	del sam
	gc.collect()

	return (tol_host, tol_rootspec, tol_mapped)

def readMapping(reads, dbs, threads, add_options, seed_length, min_score, samfile, logfile, readmapper):
	"""
	Mapping reads to reference index(es)
	"""
	input_file = " ".join(reads)
	add_options = "" if not add_options else add_options
	sam_list = []
	num_input_read = 0

	if not os.path.exists( "%s/raw_sam"%argvs.tempdir ):
		os.makedirs( "%s/raw_sam"%argvs.tempdir )

	for db in dbs:
		print_message( "Mapping to %s..." % db, argvs.silent, begin_t, logfile )

		sam_outpath = "%s/raw_sam/%s.sam"%(argvs.tempdir, os.path.basename(db))

		bash_cmd   = "set -o pipefail; set -x;"
		alnOnly_cmd= "gawk -F\\\\t '!/^@/ { print }'"
		
		# TODO: Using minimap2 on short reads requires --secondary=yes option. Therefore,
		#       the non-primary alignments (0x100) shouldn't be filtered. Also, the SAM 
		#       output of minimap2 isn't ordered by read name. Need to rewrite read-counting
		#       codes. 

		# filtered not primary alignment (0x100)
		count_cmd  = "gawk -F\\\\t '!and($2,256) && !and($2,2048) { print } END { print NR > \"%s.count\" }'"%sam_outpath
		filter_cmd = "gawk -F\\\\t '!and($2,4) { print }'"
		
		if readmapper == "minimap2":
			mapper_cmd = f"minimap2 -aL -t {threads} {add_options} {db} {input_file}"
		else:
			mapper_cmd = f"bwa mem -k{seed_length} -T{min_score} {add_options} -t{threads} {db} {input_file}"

		cmd = f"{bash_cmd} {mapper_cmd} 2>> {logfile} | {alnOnly_cmd} | {count_cmd} | {filter_cmd} > {sam_outpath}"

		if argvs.verbose: print_message( "[INFO] CMD: %s"%cmd, argvs.silent, begin_t, logfile )

		proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
		outs, errs = proc.communicate()
		exitcode = proc.poll()

		if exitcode!=0:
			print_message( "[ERROR] error occurred while running read mapping (code: %s, message: %s)."%(exitcode, errs), argvs.silent, begin_t, logfile, True )
		else:
			sam_list.append(sam_outpath)

		# get the total number of input reads
		if not num_input_read:
			with open(sam_outpath+".count") as f:
				line = f.readline()
				num_input_read = int( line.strip() )
				f.close()		
			if not num_input_read > 0:
				print_message( "[ERROR] No input reads.", argvs.silent, begin_t, logfile, True )

	#sorting sam list from small to large
	sam_list.sort(key=lambda filename: os.path.getsize(filename))

	with open(samfile, "w") as s:
		print_message( "Done mapping reads to the database(s).", argvs.silent, begin_t, logfile )
		print_message( "Merging SAM files...", argvs.silent, begin_t, logfile )
		(tol_host, tol_rootspec, tol_mapped) = mergingSAM( sam_list, s, min_score, "H", argvs.extraNM )
	s.close()

	return (tol_host, tol_rootspec, tol_mapped, num_input_read)

def print_message(msg, silent, start, logfile, errorout=0):
	message = "[%s] %s\n" % (timeSpend(start), msg)
	#loging
	with open( logfile, "a" ) as f:
		f.write( message )
		f.close()
	
	if errorout:
		sys.exit( message )
	elif not silent:
		sys.stderr.write( message )

def timeSpend( start ):
	"""
	Output time stamps.
	"""
	done = time.time()
	elapsed = done - start
	return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def getStrainSize(dbs):
	"""
	Get genome size from ann files
	
	ARGUMENT:
		dbs    LIST  input BWA indexes
	RETURN
		strain DICT  strain size in DICT with tax id as keys
		host   DICT  dict of host's taxids
	"""
	strain={}
	host={}
	for db in dbs:
		if os.path.isfile(db.rstrip(".mmi")+".ann"):
			with open(db.rstrip(".mmi")+".ann") as f:
				for line in f:
					if '|' in line:
						tmp = line.split(' ')
						r = tmp[1].split('|')
						tid = r[2]
						if 'H' in r[3]:
							host[tid] = True
						else:
							if tid in strain:
								strain[tid] += int(r[1])
							else:
								strain[tid] = int(r[1])
		# WARNING: This is not working!!!
		# The mappy library seems to have issue not reporting all seq_names.
		elif os.path.isfile(db):
			import mappy
			mmi = mappy.Aligner(fn_idx_in = db)
			for header in mmi.seq_names:
				r = header.split('|')
				tid = r[2]
				if 'H' in r[3]:
					host[tid] = True
				else:
					if tid in strain:
						strain[tid] += int(r[1])
					else:
						strain[tid] = int(r[1])
			
			if len(strain) == 0:
				print_message( "[ERROR] Not able to retrieve genome sizes from *.ann or *.mmi.", argvs.silent, begin_t, logfile, True )

	return (strain, host)

def loadBgMask(bgfiles):
	"""
	Loading background mask files
	"""
	import gzip

	bg_mask = {}
	for bgfile in bgfiles:
		try:
			with gzip.GzipFile(bgfile, 'r') as f:
				mask = json.loads(f.read())
				for ref in mask:
					maskint = int(mask[ref][2:], 16)
					if ref in bg_mask:
						bg_mask[ref] |= maskint
					else:
						bg_mask[ref] = maskint
				if argvs.verbose: print_message( "[INFO] JSON file %s loaded."%bgfile, argvs.silent, begin_t, logfile )
		except IOError:
			sys.stderr.write( "Failed to open background mask files: %s.\n"%bgfile )

	return bg_mask

def getMaskWorker(refs):
	"""
	Worker function for parsing depth file to binary mask for each reference sequences.
	"""

	mask={}
	for ref in refs:
		depfile = "%s/merged_sam/%s.sorted.depth"%(argvs.tempdir, ref)
		if os.path.isfile(depfile):
			(acc,slen,tid,tag) = ref.split('|')
			df = pd.read_csv(depfile,
				sep='\t',
				names=['ref','pos','dep'],
				usecols=['pos','dep'],
				dtype={'pos':int,'dep':int},
				index_col=['pos']
			)
			if len(df)==0: return mask
			start = int(df.head(1).index[0])
			slen = int(slen)
			df = df.reindex(range(start, slen+1), fill_value=0)
			bitmask = "".join( df['dep'].apply(lambda x: "1" if x>0 else "0").tolist() )
			mask[ref] = hex(int(bitmask,2))
	return mask

def remove_unqualified_refs(refs, res_rollup, mb, mr, ml, mc):
	qualified_refs=[]
	for ref in refs:
		(acc,leng,tid,tags) = ref.split('|')
		if mb > res_rollup[tid]["rsMRNb"] or mr > res_rollup[tid]["MR"] or ml > res_rollup[tid]["LL"] or mc > res_rollup[tid]["LC"]:
			continue
		else:
			qualified_refs.append(ref)
	
	return qualified_refs

def buildMask(refs, numthreads):
	"""
	Building bit mask for each reference sequences in parallel.
	"""

	pool = Pool(processes=numthreads)
	jobs = []
	results = []
	mask = {}

	n = 10
	chunks = [refs[i:i + n] for i in range(0, len(refs), n)]
	for chunk in chunks:
		jobs.append( pool.apply_async(getMaskWorker, (chunk,)) )

	tol_jobs = len(jobs)
	if argvs.verbose: print_message( "[INFO] Progress: %s jobs pooled."%tol_jobs, argvs.silent, begin_t, logfile )

	cnt=0
	for job in jobs:
		refmasks = job.get()
		if refmasks:
			mask.update(refmasks)
		cnt+=1
		if argvs.verbose: print_message( "[INFO] Progress: %s/%s (%.1f%%) job finished."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

	#clean up
	pool.close()

	return mask

if __name__ == '__main__':
	argvs    = parse_params( __version__ )
	begin_t  = time.time()
	sam_fp   = argvs.sam[0] if argvs.sam else ""
	samfile  = "%s/%s.pangia.sam" % (argvs.outdir, argvs.prefix) if not argvs.sam else sam_fp.name
	logfile  = "%s/%s.pangia.log" % (argvs.outdir, argvs.prefix)
	saveJson = "%s/%s.pangia.json.gz" % (argvs.outdir, argvs.prefix) if argvs.saveBg else ""
	lines_per_process = 15000
	patho_meta = {}
	genome_size = {}
	host_genome = {}
	bg_mask = {}

	#create output directory if not exists
	if not os.path.exists(argvs.outdir):
		os.makedirs(argvs.outdir)

	print_message( "Starting PanGIA %s" % __version__, argvs.silent, begin_t, logfile )

	#dependency check
	if sys.version_info < (3,0):
		print_message( "[ERROR] Python 3.0 or above is required.", argvs.silent, begin_t, logfile, True )
	if not dependency_check("bwa"):
		print_message( "[ERROR] Executable bwa not found.", argvs.silent, begin_t, logfile, True )
	if not dependency_check("gawk"):
		print_message( "[ERROR] Executable gawk not found.", argvs.silent, begin_t, logfile, True )
	if not dependency_check("samtools"):
		print_message( "[ERROR] Executable samtools not found.", argvs.silent, begin_t, logfile, True )
	if not dependency_check("parallel"):
		print_message( "[ERROR] Executable parallel not found.", argvs.silent, begin_t, logfile, True )
	if not dependency_check("minimap2"):
		print_message( "[ERROR] Executable minimap2 not found.", argvs.silent, begin_t, logfile, True )

	#prepare output object
	out_fp = sys.stdout
	outfile = "STDOUT"
	if not argvs.stdout:
		#create output directory if not exists
		if not os.path.exists(argvs.outdir):
			os.makedirs(argvs.outdir)
		ext = "fastq" if argvs.mode == "extract" else "tsv"
		tg_taxid = argvs.taxonomy if argvs.taxonomy else ""
		outfile = "%s/%s.%s%s.%s" % ( argvs.outdir, argvs.prefix, argvs.mode, tg_taxid, ext)
		out_fp = open( outfile, encoding='utf-8', mode='w')
	
	#prepare Temporary directory
	if os.path.exists(argvs.tempdir):
		print_message( "Temporary directory '%s' found. Deleting directory..."%argvs.tempdir, argvs.silent, begin_t, logfile )
		subprocess.run( "find %s -type f -exec rm {} +"%argvs.tempdir, shell=True, check=True )
		subprocess.run( "rm -rf %s"%argvs.tempdir, shell=True, check=True)
	
	os.makedirs(argvs.tempdir)

	# "cite" GNU parallel for user since it's already cited in PanGIA publication
	if "HOME" in os.environ and not os.path.exists(os.environ["HOME"]+"/.parallel/will-cite"):
		subprocess.run("mkdir -p $HOME/.parallel; touch $HOME/.parallel/will-cite", shell=True, check=True)

	# arguments info
	print_message( "Arguments and dependencies checked:", argvs.silent, begin_t, logfile )
	print_message( "    Input reads       : %s" % argvs.input,       argvs.silent, begin_t, logfile )
	print_message( "    Input SAM file    : %s" % samfile,           argvs.silent, begin_t, logfile )
	print_message( "    Input background  : %s" % (argvs.loadBg if argvs.loadBg else "None"), argvs.silent, begin_t, logfile )
	print_message( "    Save background   : %s" % (saveJson if argvs.saveBg else "None"),     argvs.silent, begin_t, logfile )	
	print_message( "    Scoring method    : %s" % argvs.scoreMethod, argvs.silent, begin_t, logfile )
	print_message( "    Scoring parameter : %s" % argvs.scoreParam,  argvs.silent, begin_t, logfile )
	print_message( "    Database          : %s" % argvs.database,    argvs.silent, begin_t, logfile )
	print_message( "    Abundance         : %s" % argvs.relAbu,      argvs.silent, begin_t, logfile )
	print_message( "    Output path       : %s" % argvs.outdir,      argvs.silent, begin_t, logfile )
	print_message( "    Prefix            : %s" % argvs.prefix,      argvs.silent, begin_t, logfile )
	print_message( "    Mode              : %s" % argvs.mode,        argvs.silent, begin_t, logfile )
	print_message( "    Specific taxid    : %s" % argvs.taxonomy,    argvs.silent, begin_t, logfile )
	print_message( "    Threads           : %d" % argvs.threads,     argvs.silent, begin_t, logfile )
	print_message( "    First #refs in XA : %s" % argvs.procAltRefs, argvs.silent, begin_t, logfile )
	print_message( "    Extra NM in XA    : %s" % argvs.extraNM,     argvs.silent, begin_t, logfile )
	print_message( "    Minimal score     : %s" % argvs.minScore,    argvs.silent, begin_t, logfile )
	print_message( "    Minimal RSNB      : %s" % argvs.minRsnb,     argvs.silent, begin_t, logfile )
	print_message( "    Minimal reads     : %s" % argvs.minReads,    argvs.silent, begin_t, logfile )
	print_message( "    Minimal linear len: %s" % argvs.minLen,      argvs.silent, begin_t, logfile )
	print_message( "    Minimal genome cov: %s" % argvs.minCov,      argvs.silent, begin_t, logfile )
	print_message( "    Minimal depth (DC): %s" % argvs.minDc,       argvs.silent, begin_t, logfile )
	print_message( "    Minimal RSDCnr    : %s" % argvs.minRsdcnr,     argvs.silent, begin_t, logfile )
	print_message( "    Aligner option    : %s" % argvs.addOptions,  argvs.silent, begin_t, logfile )
	print_message( "    Aligner seed len  : %s" % argvs.alignSeedLength,   argvs.silent, begin_t, logfile )
	print_message( "    Aligner min score : %s" % argvs.alignMinScore,     argvs.silent, begin_t, logfile )
	print_message( "    Aligner path      : %s" % dependency_check(argvs.readmapper), argvs.silent, begin_t, logfile )
	print_message( "    Samtools path     : %s" % dependency_check("samtools"), argvs.silent, begin_t, logfile )

	if argvs.debug:
		import numpy as np
		import platform
		print_message( "Environment:", argvs.silent, begin_t, logfile )
		print_message( "    Platform version  : %s" % platform.version(), argvs.silent, begin_t, logfile )
		print_message( "    Python version    : %s" % platform.python_version(), argvs.silent, begin_t, logfile )
		print_message( "    Python compiler   : %s" % platform.python_compiler(), argvs.silent, begin_t, logfile )
		print_message( "    Numpy version     : %s" % np.__version__, argvs.silent, begin_t, logfile )
		print_message( "    Pandas version    : %s" % pd.__version__, argvs.silent, begin_t, logfile )

	#load taxonomy
	print_message( "Loading taxonomy information...", argvs.silent, begin_t, logfile )
	t.loadTaxonomy( argvs.dbPath )
	print_message( "Done.", argvs.silent, begin_t, logfile )

	#load pathogen
	print_message( "Loading pathogen information...", argvs.silent, begin_t, logfile )
	patho_meta = p.loadPathogen( argvs.dbPath + "/pathogen.tsv" )
	print_message( "Done. %s pathogens loaded."%len(patho_meta), argvs.silent, begin_t, logfile )

	#load uniqueness
	print_message( "Loading taxonomic uniqueness information...", argvs.silent, begin_t, logfile )
	uniq_meta = s.loadUniqueness( argvs.dbPath + "/uniqueness.tsv" )
	print_message( "Done. %s taxonomic uniqueness loaded."%len(uniq_meta), argvs.silent, begin_t, logfile )

	#load strain genome length
	print_message( "Loading sizes of genomes...", argvs.silent, begin_t, logfile )
	if argvs.database:
		genome_size, host_genome = getStrainSize( argvs.database )
		print_message( "Done. %s target and %s host genome(s) loaded."%(len(genome_size),len(host_genome)), argvs.silent, begin_t, logfile )

	# if reads provided
	if argvs.input:
		print_message( "Running read-mapping...", argvs.silent, begin_t, logfile )
		(tol_host, tol_rootspec, tol_mapped, tol_reads) = readMapping( argvs.input, argvs.database, argvs.threads, argvs.addOptions, argvs.alignSeedLength, argvs.alignMinScore, samfile, logfile, argvs.readmapper)
		print_message( "Logfile saved to %s." % logfile, argvs.silent, begin_t, logfile )
		print_message( "Done. Mapped SAM file saved to %s." % samfile, argvs.silent, begin_t, logfile )
		print_message( "Total number of input reads: %s"%tol_reads, argvs.silent, begin_t, logfile )
		print_message( "Total number of mapped reads: %s"%tol_mapped, argvs.silent, begin_t, logfile )

		if not (tol_mapped-tol_host-tol_rootspec):
			print_message( "No target reads found.", argvs.silent, begin_t, logfile )
			#prepare Temporary directory
			if os.path.exists(argvs.tempdir) and not argvs.verbose and not argvs.keepTemp:
				print_message( "Cleaning up temporary directory '%s'..."%argvs.tempdir, argvs.silent, begin_t, logfile )
				subprocess.run( "find %s -type f -exec rm {} +"%argvs.tempdir, shell=True, check=True )
				subprocess.run( "rm -rf %s"%argvs.tempdir, shell=True, check=True )
				print_message( "Done.", argvs.silent, begin_t, logfile )
			print_message( "PanGIA stopped.", argvs.silent, begin_t, logfile, True )

		print_message( "Total number of host reads: %s (%.2f%%)" %(tol_host, tol_host/tol_mapped*100), argvs.silent, begin_t, logfile )
		print_message( "Total number of ignored reads (cross superkingdom): %s (%.2f%%)" %(tol_rootspec, tol_rootspec/tol_mapped*100), argvs.silent, begin_t, logfile )
		sam_fp = open( samfile, "r" )

	if argvs.mode == 'class':
		print_message( "Classifying reads... ", argvs.silent, begin_t, logfile ) 
		processSAMfileReadClass( sam_fp, out_fp, argvs.taxonomy )
		print_message( "Done classifying reads. Results printed to %s." % outfile, argvs.silent, begin_t, logfile )

	elif argvs.mode == 'extract':
		print_message( "Extracting reads... ", argvs.silent, begin_t, logfile )
		processSAMfileReadExtract( os.path.abspath(samfile), out_fp, argvs.taxonomy, argvs.threads )
		print_message( "Done extracting reads to %s." % outfile, argvs.silent, begin_t, logfile )

	else:
		print_message( "Processing SAM file... ", argvs.silent, begin_t, logfile )
		(res, mapped_r_cnt) = processSAMfile( os.path.abspath(samfile), argvs.threads, lines_per_process)
		print_message( "Done processing SAM file, %s alignment(s)."% (mapped_r_cnt), argvs.silent, begin_t, logfile )

		print_message( "Rolling up taxonomies...", argvs.silent, begin_t, logfile )
		res_rollup = taxonomyRollUp(res, patho_meta, mapped_r_cnt, argvs.minRsnb, argvs.minReads, argvs.minLen, argvs.minCov, argvs.minDc)
		print_message( "Done.", argvs.silent, begin_t, logfile )

		if argvs.scoreMethod != "standalone":
			#load background
			if argvs.loadBg:
				print_message( "Loading background information...", argvs.silent, begin_t, logfile )
				bg_mask = loadBgMask(argvs.loadBg)
				print_message( "Done. %s background taxonomies loaded."%len(bg_mask), argvs.silent, begin_t, logfile )

			if bg_mask:
				print_message( "Calculating score against background...", argvs.silent, begin_t, logfile )
				refs = list(res.keys())
				res_rollup = s.scoreBg(res_rollup, refs, bg_mask, argvs.tempdir, argvs.threads, argvs.minRsnb, argvs.minReads, argvs.minLen, argvs.minCov, argvs.verbose, argvs.debug)
				print_message( "Done calculating.", argvs.silent, begin_t, logfile )

		if argvs.scoreMethod != "bg":
			print_message( "Calculating standalone score...", argvs.silent, begin_t, logfile )
			res_rollup = scoreStandalone(res_rollup, uniq_meta, argvs.scoreParam)
			print_message( "Done calculating.", argvs.silent, begin_t, logfile )

		# save res_rollup to JSON
		if argvs.saveBg:
			print_message( "Calculating bitmasks for mapped references...", argvs.silent, begin_t, logfile )
			with gzip.GzipFile(saveJson, 'w') as f:
				# get all possible refs
				refs = list(res.keys())
				# ignored if the taxa of refs do not reported finally
				#refs = remove_unqualified_refs(refs, res_rollup, argvs.minRsnb, argvs.minReads, argvs.minLen, argvs.minCov)
				# create bitmasks
				bitmask = buildMask(refs, argvs.threads)
				print_message( "Done.", argvs.silent, begin_t, logfile )
				print_message( "Saving bitmasks in JSON format...", argvs.silent, begin_t, logfile )
				f.write(json.dumps(bitmask, indent=4).encode('utf-8'))
			print_message( "Done.", argvs.silent, begin_t, logfile )

		if argvs.mode == 'report':
			outputResultsAsReport( res_rollup, out_fp, argvs.relAbu, argvs.taxonomy, argvs.reportFields, argvs.scoreMethod, argvs.minScore, argvs.minCov, argvs.minRsnb, argvs.minReads, argvs.minLen, argvs.minDc, argvs.minRsdcnr, argvs.displayAll )
		elif argvs.mode == 'lineage':
			outputResultsAsLineage( res_rollup, out_fp, argvs.relAbu, argvs.mode, argvs.scoreMethod, argvs.minScore, argvs.minRsnb, argvs.minReads, argvs.minLen, argvs.minDc, argvs.displayAll )

		print_message( "Done taxonomy profiling; %s results saved to %s." % (argvs.mode, outfile), argvs.silent, begin_t, logfile )

	#prepare Temporary directory
	if os.path.exists(argvs.tempdir) and not argvs.verbose and not argvs.keepTemp:
		print_message( "Cleaning up temporary directory '%s'..."%argvs.tempdir, argvs.silent, begin_t, logfile )
		subprocess.run( "find %s -type f -exec rm {} +"%argvs.tempdir, shell=True, check=True )
		subprocess.run( "rm -rf %s"%argvs.tempdir, shell=True, check=True )
		print_message( "Done.", argvs.silent, begin_t, logfile )