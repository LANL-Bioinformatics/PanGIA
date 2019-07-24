#!/usr/bin/env python
import argparse as ap, textwrap as tw
import taxonomy as t
import gzip
import sqlite3
import sys
import os
import subprocess
from re import search

def parse_params():
	p = ap.ArgumentParser( prog='pangia_db_build.py', description="""Convert sequences to PanGIA v2 format.""" )

	p.add_argument(
		'-i','--input', metavar='[FILE(S)]', type=str, nargs='*',
		help="input one or multiple FASTA file(s) (default: none)")

	p.add_argument(
		'-a', '--accPrefixOnly', metavar='[STR]', type=str, default=None,
	    help="""Parsed sequences with particular prefix in the accession number only [default: None]""")

	p.add_argument(
		'-dp', '--dbPath', metavar='[PATH]', type=str, default=None,
	    help="""Path of taxonomy databases [default: <BIN_DIR>/database/]""")

	p.add_argument(
		'-sdb','--sqlitedb', metavar='[FILE]', nargs='?',
		help="Custom SQLite3 taxonomy db [default: <PREFIX>.custom.db]")

	p.add_argument(
		'-asr', '--assemblySummaryRefseq', metavar='[PATH]', type=str, default=None,
	    help="""Use assembly_summary_refseq.txt file to define strain [default: None]""")

	p.add_argument(
		'-asd', '--assemblyRefseqDir', metavar='[PATH]', type=str,
	    help="""local directory for FASTA files [default: <PREFIX>_refseq_genomes]""")

	p.add_argument(
		'-rc', '--refseqCategoryOnly', nargs='*', metavar='[STR]', type=str,
	    help="""Only parse specified refseq category(ies), like 'reference' or 'representative'. Only works with '--assemblySummaryRefseq' option. [default: None]""")

	p.add_argument(
		'-al', '--assemblyLevelOnly', nargs='*', metavar='[STR]', type=str,
	    help="""Only parse specified assembly level(s), like 'complete'. Only works with '--assemblySummaryRefseq' option. [default: None]""")

	p.add_argument(
		'-t','--tag', metavar='[STR]', type=str, 
		help="Specify a TAG and added to all input sequences. Use 'H' when you specify host reference sequences. The initial of superkingdom name will be used by default, for example 'B' for Bacteria and 'V' for Viruses.")

	p.add_argument(
		'-ssf', '--strainSeparateFasta', action="store_true",
	    help="""Each FASTA file represents a strain. This option is ALWAYS TRUE if --assemblySummaryRefseq is used. [default: FALSE]""")

	p.add_argument(
		'-sr', '--splitDatabaseByRank', metavar='[RANK]', type=str, default="superkingdom",
	    help="""Split database by rank [default: superkingdom]""")

	p.add_argument(
		'-x', '--taxid', metavar='[TAXID]', nargs='*', type=str, default=None,
	    help="""Only process sequences belong to specified taxid(s) [default: None]""")

	p.add_argument(
		'-s','--skipSeqType', metavar='[STR]', nargs='*', type=str,
		help="""Skip specific type(s) of sequence: 'p' for plasmid and 'g' for phage. [default: None]""")

	p.add_argument(
		'-ngl','--notGottchaDbList', action="store_true",
		help="""Not generate gottcha_db input list. [default: False]""")

	p.add_argument(
		'-p', '--prefix', metavar='[STR]', type=str, default="pangia",
	    help="""Output file [default: pangia]""")

	p.add_argument(
		'-rs', '--redundantStrain', action="store_true",
	    help="""Allow redundant strains [default: FALSE]""")

	p.add_argument(
		'-iso', '--includeIsolate', action="store_true",
	    help="""Different isolates will be treated as differnet strains [default: FALSE]""")

	args_parsed = p.parse_args()

	if not args_parsed.dbPath:
		bin_dir = os.path.dirname(os.path.realpath(__file__))
		args_parsed.dbPath = bin_dir + "/database"

	if not args_parsed.assemblyRefseqDir:
		args_parsed.assemblyRefseqDir = args_parsed.prefix + "_refseq_genomes"

	if not args_parsed.sqlitedb:
		args_parsed.sqlitedb = args_parsed.prefix + ".custom.db"

	if args_parsed.assemblySummaryRefseq:
		args_parsed.strainSeparateFasta = True

	if not args_parsed.input and not args_parsed.assemblySummaryRefseq:
		sys.stderr.write( "[ERROR] Neither --input nor --assemblySummaryRefseq wer specified.\n" )
		sys.exit(1)

	return args_parsed

def getAccFromHeader( full_header ):
	acc = 'unknown'
	header = full_header.split(" ")[0]

	# old RefSeq header:
	# >gi|546540813|ref|NZ_CBML010000001.1| Clostridium chauvoei JF4335 WGS project CBML000000000 data, contig 00016
	if 'ref|' in header:
		m = search( 'ref\|([^\|]+)', header)
		try:
			acc = m.group(1)
		except:
			sys.stderr.write( "\n[WARNING] seems 'ref', but failed to parse; header: %s\n" % header )
	elif 'gb|' in header:
		m = search( 'gb\|([^\|]+)', header)
		try:
			acc = m.group(1)
		except:
			sys.stderr.write( "\n[WARNING] seems 'gb', but failed to parse; header: %s\n" % header )
	# >CBML010
	elif '|' not in header:
		acc = header[1:]
	else:
		sys.stderr.write( "\n[WARNING] Unrecognized format; header: %s\n" % header )

	return acc

def generateCustomStrain(taxid, name, desc, c, asr_strain, asr_isolate, include_iso):
	#get custom strain name
	(str_name, iso_name) = getStrainName(taxid, name, desc, asr_strain, asr_isolate)

	# get full custom strain name
	cus_name = name
	if str_name and not str_name in cus_name:
		cus_name += " strain " + str_name

	if '[strain|' in desc:
		m = search( '\[strain\|([^\]]+)\]', desc )
		if m and not m.group(1) in name:
			cus_name = m.group(1)

	if iso_name and include_iso and not iso_name in cus_name:
		cus_name += " isolate " + iso_name

	#check custom name's existence
	c.execute( 'SELECT * FROM taxonomy_custom WHERE P_TAXID="%s" and NAME="%s"' % (taxid, cus_name) )
	tax_rec = c.fetchone()

	if tax_rec:
		cus_taxid = tax_rec[0]
	else:
		cus_taxid = getNewCusTaxonomy( c, taxid, cus_name )

	return cus_taxid, cus_name, str_name, iso_name

def getNewCusTaxonomy( c, taxid, cus_name, retry=0 ):
	#insert new custom name
	cus_taxid = ""
	str_id = ""
	depth = ""
	c.execute('SELECT * FROM taxonomy_custom WHERE P_TAXID="%s" ORDER BY STR_ID DESC' % taxid )
	tax_rec = c.fetchone()
	# if parent taxid (usually a species) already in the database
	if tax_rec:
		str_id = tax_rec[2]+1
		depth  = tax_rec[4]
	else:
		str_id = 1
		depth  = int(t.taxid2depth(taxid)) + 1
	
	cus_taxid = "%s.%s" % (taxid, str_id)
	c.execute('INSERT INTO taxonomy_custom VALUES ("%s", "%s", %s, "%s", %s, "%s")' % ( cus_taxid, taxid, str_id, cus_name, depth, "no rank") )

	return cus_taxid

def checkNewTaxid( c, tid, ptid, name, str_name="", iso_name="", aacc="", filename="" ):
	"""
	ARGUMENTS
	    c
	    tid
	    name
	    str_name
	    iso_name
	    aacc
	    filename
	RETURN
		BOOL
	"""

	c.execute('SELECT * FROM strain_taxonomy WHERE TAXID="%s"' % tid )
	tax_rec = c.fetchone()
	# if another assembly is named the same strain, skip this file (assembly)
	if tax_rec:
		sys.stderr.write( "[INFO] Strain exists (Taxid: %s, Name: %s).\n" % (tid, tax_rec[2]) )
		return False
	# else, add this strain to the taxanomy database
	else:
		#    TAXID      CHAR(30)    NOT NULL, /* taxid */
		#    P_TAXID    CHAR(30)    NOT NULL, /* parent taxid */
		#    NAME       CHAR(200)   NOT NULL, /* taxonomy name */
		#    STR_NAME   CHAR(200)   NULL,     /* strain name */
		#    ISO_NAME   CHAR(200)   NULL,     /* isolate name */
		#    ASSEM_ACC  CHAR(50)    NULL,     /* assembly acc# */
		#    FILE       CHAR(50)    NULL,     /* file name */
		#    PRIMARY KEY (TAXID)
		c.execute('INSERT INTO strain_taxonomy VALUES ("%s", "%s", "%s", "%s", "%s", "%s", "%s")' % ( 
			tid, 
			ptid, 
			name, 
			str_name, 
			iso_name, 
			aacc, 
			filename
			)
		)
		return True

def getStrainName(taxid, name, desc, asr_strain, asr_isolate):
	"""
	ARGUMENTS
	    taxid
		name
		desc
		asr_strain
		asr_isolate
	RETURN
		str_name
		iso_name
	"""
	str_name = ""
	iso_name = ""

	# get strain info from assembly summary report
	if asr_strain:
		str_name = asr_strain
		iso_name = asr_isolate

	#replace ", " with space
	mod_desc = desc.replace( ", ", "  " )
	mod_desc = desc.replace( '"', "'" )

	#if custom strain name exists, use it; otherwise, try to identify it from desc
	if '[isolate|' in desc:
		m = search( '\[isolate\|([^\]]+)\]', mod_desc )
		if m and not m.group(1) in name: iso_name = m.group(1)
	else:
		if 'isolate' in desc:
			m = search( '(isolate \S+)', mod_desc )
			if m and not m.group(1) in name: iso_name = m.group(1)
		elif 'clone' in desc:
			m = search( '(clone \S+)', mod_desc )
			if m and not m.group(1) in name: iso_name = m.group(1)
	
	#for custom isolate name exists, use it; otherwise, try to identify it from desc
	if '[strain|' in desc:
		m = search( '\[strain\|([^\]]+)\]', mod_desc )
		if m and not m.group(1) in name: str_name = m.group(1)
	else:
		if 'strain' in desc:
			m = search( '(strain \S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)
		elif 'strain:' in desc:
			m = search( '(strain: \S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)
		elif ' str ' in desc:
			m = search( ' str (\S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)
		elif ' str. ' in desc:
			m = search( ' str\. (\S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)
		elif 'sp.' in desc:
			m = search( '(sp\. \S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)
		elif '/' in desc: # for viruses especially
			m = search( '(\S+\/\S+\/\S+\/\S+)', mod_desc )
			if m and not m.group(1) in name: str_name += " "+m.group(1)

		str_name = str_name.rstrip(',')

		#more strain name search if all of above failed
		if not str_name:
			if "ATCC " in desc:
				m = search( '(ATCC \S+)', mod_desc )
				if m and not m.group(1) in name: str_name += " "+m.group(1)
			if "NCTC " in desc:
				m = search( '(NCTC \S+)', mod_desc )
				if m and not m.group(1) in name: str_name += " "+m.group(1)
			elif "genome assembly " in desc:
				m = search( 'genome assembly (\S+)', mod_desc )
				if m and not m.group(1) in name: str_name += " "+m.group(1)	
			elif "," in desc:
				if name in desc:
					m = search( '%s (\S+),'%name, desc )
					if m and not m.group(1) in name: str_name += " "+m.group(1)
				elif name in desc:
					m = search( '%s (\S+)'%name, desc )
					if m and not m.group(1) in name: str_name += " "+m.group(1)
				elif name in desc:
					m = search( '%s ([^,]+),'%name, desc )
					if m and not m.group(1) in name: str_name += " "+m.group(1)
				else:
					m = search( '(\S+),', desc )
					if m and not m.group(1) in name:
						if not m.group(1).startswith('('):
							str_name += " "+m.group(1).rstrip(')')
						else:
							str_name += " "+m.group(1)
		
		# desperate mode
		if not str_name:
			str_name = mod_desc

	return (str_name, iso_name)

def initTaxaDB(c):
	c.execute("""
	CREATE TABLE IF NOT EXISTS taxonomy_custom(
	   TAXID      CHAR(30)    NOT NULL,
	   P_TAXID    CHAR(20)    NOT NULL,
	   STR_ID     INT         NOT NULL,
	   NAME       CHAR(200)   NOT NULL,
	   DEPTH      INT         NOT NULL,
	   RANK       CHAR(20)    NOT NULL,
	   PRIMARY KEY (TAXID)
	);
	""")

	c.execute("""
	CREATE TABLE IF NOT EXISTS strain_taxonomy(
	   TAXID      CHAR(30)    NOT NULL, /* taxid */
	   P_TAXID    CHAR(30)    NOT NULL, /* parent taxid */
	   NAME       CHAR(200)   NOT NULL, /* taxonomy name */
	   STR_NAME   CHAR(200)   NULL,     /* strain name */
	   ISO_NAME   CHAR(200)   NULL,     /* isolate name */
	   ASSEM_ACC  CHAR(50)    NULL,     /* assembly acc# */
	   FILE       CHAR(50)    NULL,     /* file name */
	   PRIMARY KEY (TAXID)
	);
	""")

	conn.commit()

def processASR( c, tid2lineage, asr, rc, al, asr_dir, target_tid, tag, flag_red_str, flag_ind_iso=0 ):
	"""
	ARGVS:
	   c             OBJ   sqlite3 connection obj
	   asr           STR   asr file
	   rc            LIST  refseq category
	   al            LIST  assembly level
	   asr_dir       STR   asr local directory
	   target_tid    LIST  target taxid
	   tag           STR   tag
	   flag_red_str  BOOL  flag for allowing redundant strains
	   flag_ind_iso  BOOL  flag for including isolate in strain name
	RETURN:
	   asr_info      DICT
	"""

	# init vars
	asr_info = t._autoVivification()
	lineage  = t._autoVivification()
	cnt_q    = 0 #qualified assembly
	cnt_tol  = 0 #total genomes

	# parsing assembly_summary_refseq.txt file
	with open(asr) as f:
		for line in f:
			if line.startswith('#'):
				continue
			else:
				cnt_tol += 1

			line = line.strip('\n')
			
			#split each line in assembly_summary_refseq.txt:
			#  0- 4  assembly_accession    bioproject      biosample         wgs_master           refseq_category
			#  5- 9  taxid                 species_taxid   organism_name     infraspecific_name   isolate 
			# 10-14  version_status        assembly_level  release_type      genome_rep           seq_rel_date
			# 15-19  asm_name              submitter       gbrs_paired_asm   paired_asm_comp      ftp_path        
			# 20-21  excluded_from_refseq  relation_to_type_material
			tmp = line.split('\t')
			filename = ""
			local_path = ""

			sys.stderr.write( "[INFO] Processing: %s..."%tmp[0] )

			# try to get taxonomy lineage
			try:
				lineage = t.taxid2lineageDICT(tmp[5], 1, 1)
			except:
				sys.stderr.write( "skipped. Removed TaxID (%s) found for %s.\n" % (tmp[5], tmp[0]) )
				continue
			# SKIPPING following records
			# 1) not specified refseq_category, example: reference 
			if rc:
				flag=0
				for cate in rc:
					if cate.lower() in tmp[4].lower():
						flag=1
						break
				if not flag:
					sys.stderr.write( "skipped. Not belongs to specific refseq_category %s.\n"%rc )
					continue
			# 2) not specified assembly_level, example: complete
			if al:
				flag=0
				for l in al:
					if l.lower() in tmp[11].lower():
						flag=1
						break
				if not flag:
					sys.stderr.write( "skipped. Not belongs to specific assembly_level %s.\n"%al )
					continue
			# 3) assemblies that marked "excluded_from_refseq"
			if tmp[20]:
				sys.stderr.write( "skipped. Marked as excluded_from_refseq.\n" )
				continue
			# 4) not belongs to specified tax id
			if tmp[5] and target_tid:		
				flag_tid_in_lineage=0
				for t_tid in target_tid:
					if not flag_tid_in_lineage:
						for t_rank in lineage:
							if lineage[t_rank]['taxid'] == t_tid:
								flag_tid_in_lineage=1
								break
				if not flag_tid_in_lineage:
					sys.stderr.write( "skipped. Not belongs to specified tid.\n" )
					continue
			# 5) sequence location is not available
			if tmp[19] != "na":
				if tmp[19].startswith('ftp'):
					# real filename of assembly
					filename = tmp[19].split('/')[-1] + "_genomic.fna.gz"
					# directory structure is retained to keep local filesystem healthy
					local_path = asr_dir + tmp[19].split('genomes')[1]
				else:
					filename = tmp[19].split('/')[-1]
					local_path = "/".join(tmp[19].split('/')[:-1])
			else:
				sys.stderr.write( "skipped. No URL provided.\n" )
				continue
			# 6) not an unique strain
			tid = tmp[5]
			rank = t.taxid2rank(tid)
			name = t.taxid2name(tid)
			if rank != "strain" and rank != "unknown":
				# generate custom strain if not exists and add strains to database
				(cus_str_taxid, cus_str_name, str_name, iso_name) = generateCustomStrain( tid, name, '', c, tmp[8].replace("strain=",""), tmp[9], flag_ind_iso )
				#add custom strain to lineage
				lineage["strain"]['taxid'] = cus_str_taxid
				lineage["strain"]['name'] = cus_str_name

			flag_new_strain = checkNewTaxid(c, lineage["strain"]['taxid'], "", lineage["strain"]['name'], tmp[8].replace("strain=",""), tmp[9], tmp[0], filename )
			if not flag_new_strain and not flag_red_str:
				sys.stderr.write( "skipped. Not an unique strain.\n" )
				continue

			sys.stderr.write( "qualified.\n" )

			# download sequence files if not available at local directory
			if not os.path.isfile( local_path+"/"+filename ):
				url = tmp[19]+"/"+filename
				wget( url, local_path )
			else:
				sys.stderr.write( "[INFO] Found local file: %s.\n"%(local_path+"/"+filename) )

			# tag
			if not tag:
				tag = lineage["superkingdom"]['name'][0]

			# use assembly_accession as key
			cnt_q += 1
			tid = lineage["strain"]['taxid']
			asr_info[tmp[0]]['taxid']         = tid
			asr_info[tmp[0]]['full_str_name'] = lineage["strain"]['name']
			asr_info[tmp[0]]['ftp_path']      = tmp[19]
			asr_info[tmp[0]]['local_path']    = local_path
			asr_info[tmp[0]]['filename']      = filename
			asr_info[tmp[0]]['type_material'] = True if tmp[21] else False
			asr_info[tmp[0]]['cate_tag']      = tag
			tid2lineage[tid] = lineage
	f.close()
	return asr_info, tid2lineage, cnt_q, cnt_tol

def processFASTA( filename, asm_asr_info, tid2lineage, tag, flag_str_sep_fa, accPrefixOnly, skipSeqType, output, flag_ind_iso=False, flag_no_gdb_list=0, flag_redundant_str=False ):
	"""
	ARGVS:
	   filename        STR   FASTA filename
	   asm_asr_info    DICT  asr_info[asm]
	   tag             STR   tag
	   flag_str_sep_fa BOOL  all sequences in the file belong to one strain
	   accPrefixOnly   STR   process sequences with accession number start with particular prefix only
	   skipSeqType     LIST  skip type(s) of sequence
	   output          OBJ   file object for output FASTA
	RETURN:
	   cnt             INT   total # of output sequences
	   cnt_tol         INT   total # of input sequences
	"""

	#init vars
	(acc, rtype, rdesc, rank, rtid, rname, rseq, rnote) = ("", "", "", "", "", "", "", "")
	ignore = False
	lineage = {}
	gdb_list = {}
	cnt_tol = 0
	cnt = 0

	# open files
	if os.path.exists(filename):
		if filename.endswith('gz'):
			infile = gzip.open(filename, mode='rt', encoding='utf-8')
		else:
			infile = open(filename, 'r')
	else:
		sys.stderr.write( "[ERROR] %s doesn't exists.\n"%filename )
		exit

	# processing FASTA files
	for line in infile:
		if line.startswith(">"):
			cnt_tol+=1
			#print last sequences

			if acc and not ignore:
				cnt+=1
				output.write( ">%s|%s|%s|%s [%s|%s] %s\n%s\n" % (acc, len(rseq), rtid, rtype, rank, rname, rdesc, rseq) )
				if rank == 'unknown': sys.stderr.write( "\n[WARNING] Unknown sequence: >%s|%s|%s|%s %s\n" % (acc, len(rseq), rtid, rtype, rdesc) )
				#reset taxa info for sequences

			if not flag_str_sep_fa:
				(acc, rtype, rdesc, rank, rtid, rname, rseq, rank) = ("", "", "", "", "", "", "", "")
			else:
				(acc, rtype, rseq) = ("", "", "")

			#parse header
			temp = line.strip().split(' ')
			header = temp[0]
			acc = getAccFromHeader( header )
			rdesc = " ".join(temp[1:])

			#somehow failed to parse acc from header
			if not acc or acc=="unknown":
				sys.stderr.write( "[WARNING] Unknown acc#: >%s|%s|%s|%s %s\n" % (acc, len(rseq), rtid, rtype, rdesc) )
				ignore = True
				continue

			# check prefix
			if accPrefixOnly:
				if not acc.startswith(accPrefixOnly):
					sys.stderr.write( "[WARNING] Ignored, acc prefix not match: >%s|%s|%s|%s %s\n" % (acc, len(rseq), rtid, rtype, rdesc) )
					ignore = True
					continue

			# get taxonomy info
			if "taxid" in asm_asr_info:
				#asm file is provided
				rtid = asm_asr_info["taxid"]
				rname = asm_asr_info["full_str_name"]
				rank = "strain"
				lineage = tid2lineage[rtid]
			else:
				if not rtid:
					# for custom acc#, parse [taxid|xxxx]
					if "[taxid|" in rdesc:
						m = search( '\[taxid\|(\d+)\]', rdesc)
						try:
							rtid = m.group(1)
						except:
							sys.stderr.write( "[WARNING] Failed to parse taxid from [taxid|xxxx].\n")
					else:
						# convert from acc2taxid
						rtid = t.acc2taxid(acc)

					if rtid and t.taxidStatus( rtid ) != "invalid":
						rank = t.taxid2rank( rtid )
						name = t.taxid2name( rtid )
					else: # skip the reference if not able to find a tax id
						ignore = True
						sys.stderr.write( "[WARNING] Ignored. Failed mapping to TaxID.\n")
						continue

					lineage = t.taxid2lineageDICT(rtid, 1, 1)

					# if rtid is not strain level taxonomy
					(str_name, iso_name) = ("","")
					if rank != "strain" and rank != "unknown":
						(new_taxid, new_name, str_name, iso_name) = generateCustomStrain( rtid, name, rdesc, c, '', '', flag_ind_iso )
						#write custom tax
						rtid  = new_taxid
						rname = new_name
						rank  = "strain"
						lineage[rank]['taxid'] = rtid
						lineage[rank]['name'] = rname
					else: # if rtid is strain level tax
						rname = name

					#check if the strain is unique
					flag_new_strain = checkNewTaxid(c, rtid, "", rname, str_name, iso_name, "", filename )

					if not rtid in tid2lineage:
						tid2lineage[rtid] = lineage

					if not flag_new_strain and not flag_redundant_str:
						sys.stderr.write( "[WARNING] Ignored, strain exists.\n")
						ignore = True
						continue

			# tags for different sequence types
			if 'cate_tag' in asm_asr_info:
				rtype = asm_asr_info['cate_tag']
			elif tag:
				rtype = tag
			else:
				rtype = lineage["superkingdom"]['name'][0]

			if 'plasmid' in line or 'Plasmid' in line:
				rtype += 'p' #plasmid
			elif 'phage' in line or 'Phage' in line:
				rtype += 'g' #Phage
			else:	
				rtype += 'c' #chromosome

			if skipSeqType and rtype in skipSeqType:
				sys.stderr.write( "[WARNING] Ignored, unwanted type.\n")
				ignore = True
				continue

			# PLASMID ONLY: if a plasmid doesn't belong to a species or up, roll it up to species
			#if 'p' in rtype and ( rank == 'strain' or rank == 'no rank' ):
			#	rank = "species"
			#	rtid = t.taxid2taxidOnRank(rtid, rank)
			#	rname = t.taxid2name(rtid)

			# hash for gotthca_db input list
			if not flag_no_gdb_list:
				gdb_list[filename] = rtid
			
			ignore = False
		else:
			if ignore: continue
			line = line.strip()
			if not line: continue #empty line
			rseq += line
	#last data
	if acc and not ignore:
		cnt+=1
		output.write( ">%s|%s|%s|%s [%s|%s] %s\n%s\n" % (acc, len(rseq), rtid, rtype, rank, rname, rdesc, rseq) )

	infile.close()
	return (cnt, cnt_tol, gdb_list)

def generateGottchaDbList(gdb_list, output_gdbl_buffer, tid2lineage):
	for fn in gdb_list:
		tid = gdb_list[fn]
		ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
		for rank in ranks:
			text = "%s\t%s\tlinear\t%s\n"%( tid2lineage[tid][rank]['name'], fn, tid )
			if not rank in output_gdbl_buffer:
				output_gdbl_buffer[rank] = text
			else:
				output_gdbl_buffer[rank] += text
	return output_gdbl_buffer

def wget(url, lpath="./"):
	subprocess.run("mkdir -p %s"%lpath, shell=True, check=True)
	subprocess.run("wget '%s' -P %s"%(url, lpath), shell=True, check=True)

if __name__ == '__main__':
	argvs    = parse_params()

	cnt         = 0
	cnt_fa      = 0 # accepted genomes
	cnt_fa_tol  = 0 # total genomes
	cnt_seq     = 0
	cnt_seq_tol = 0
	tid2lineage = t._autoVivification()
	output      = sys.stdout
	output_gdbl_buffer = {}

	# loading taxonomy
	sys.stderr.write( "Loading taxonomy..." )
	t.loadTaxonomy( argvs.dbPath )
	sys.stderr.write( "completed.\n" )

	# init the taxonomy sqlite3 db file
	conn = sqlite3.connect(argvs.sqlitedb)
	conn.isolation_level = None
	c = conn.cursor()
	initTaxaDB(c)

	# create output file
	output_fn = "%s.fasta"%argvs.prefix
	sys.stderr.write( "The script will save output sequences to %s.\n"%output_fn )
	output = open( "%s.fasta"%argvs.prefix, "w")

	sys.stderr.write( "Start processing...\n" )

	# input from assembly summary refseq report
	if argvs.assemblySummaryRefseq and os.path.isfile(argvs.assemblySummaryRefseq):
		sys.stderr.write( "Parsing assembly summary refseq report...\n" )
		(asr_info, tid2lineage, cnt_fa, cnt_fa_tol) = processASR( c, tid2lineage, argvs.assemblySummaryRefseq, argvs.refseqCategoryOnly, argvs.assemblyLevelOnly, argvs.assemblyRefseqDir, argvs.taxid, argvs.tag, argvs.redundantStrain, argvs.includeIsolate )
		sys.stderr.write( "Done. %s/%s genome(s) qualified.\n"%(cnt_fa, cnt_fa_tol) )
		
		sys.stderr.write( "Parsing genomes sequences...\n" )
		for asm in asr_info:
			input_filename = asr_info[asm]['local_path']+"/"+asr_info[asm]['filename']
			sys.stderr.write( "Parsing %s from %s...\n"%(asm,input_filename) )
			(cnt_fa_seq, cnt_fa_seq_tol, gdb_list) = processFASTA( input_filename, asr_info[asm], tid2lineage, argvs.tag, True, argvs.accPrefixOnly, argvs.skipSeqType, output, argvs.includeIsolate )
			output_gdbl_buffer = generateGottchaDbList(gdb_list, output_gdbl_buffer, tid2lineage)
			cnt_seq += cnt_fa_seq
			cnt_seq_tol += cnt_fa_seq_tol
			cnt+=1
			sys.stderr.write( "Done. %s/%s sequence(s) parsed.\n"%(cnt_fa_seq, cnt_fa_seq_tol) )
			sys.stderr.write( "%s/%s qualified genome(s) parsed.\n"%(cnt, cnt_fa) )

		if not argvs.notGottchaDbList:
			sys.stderr.write( "Generating gottcha_db input lists...\n" )
			ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
			for rank in ranks:
				with open( "%s.gottcha_db.%s.list"%(argvs.prefix, rank), "w") as f:
					f.write( output_gdbl_buffer[rank] )
				f.close()
			sys.stderr.write( "Done.\n" )

	# input from FASTA files
	if argvs.input:
		for input_filename in argvs.input:
			cnt_fa += 1
			if os.path.isfile(input_filename):
				sys.stderr.write( "Parsing %s...\n"%input_filename )
				(cnt_fa_seq, cnt_fa_seq_tol, gdb_list) = processFASTA( input_filename, {}, tid2lineage, argvs.tag, argvs.strainSeparateFasta, argvs.accPrefixOnly, argvs.skipSeqType, output, argvs.includeIsolate, argvs.notGottchaDbList, argvs.redundantStrain )
				cnt += 1
				cnt_seq += cnt_fa_seq
				cnt_seq_tol += cnt_fa_seq_tol
				sys.stderr.write( "Done. %s/%s sequence(s) parsed.\n"%(cnt_fa_seq, cnt_fa_seq_tol) )
			else:
				sys.stderr.write( "[WARNING] Skipped: file is not available: %s\n"%input_filename )

	# output custom database from sqlite3 to tsv
	with open(argvs.sqlitedb+".tsv", "w") as tsv_output:
		try:
			c.execute( 'SELECT * FROM taxonomy_custom ORDER BY NAME' )
			for row in c.fetchall():
				tsv_output.write( "%s\t%s\t%s\t%s\t%s\n"%(row[0], row[4], row[1], row[5], row[3]) )
		except sqlite3.Error as e:
			sys.stderr.write( "[ERROR] An error occurred: %s\n" % e.args[0] )

	conn.close()
	sys.stderr.write( "\nComplete! %d/%d qualified genome file(s); %d/%d sequence(s) processed. \n" % (cnt, cnt_fa, cnt_seq, cnt_seq_tol) )