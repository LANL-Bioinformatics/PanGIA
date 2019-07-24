#!/usr/bin/env python3

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = ""
__date__      = ""
__copyright__ = ""

import argparse as ap
import sys
import os
from re import search, findall
import taxonomy as t

def parse_params():
	p = ap.ArgumentParser(prog='pangia_report2lineage.py', description="""Convert PanGIA report to lineage format""")

	p.add_argument('-d', '--database',
			metavar='[BWA_INDEX]', type=str, default=None,
					help="Name/path of readmapper's index [default: None]")

	p.add_argument('-dp', '--dbPath',
			metavar='[PATH]', type=str, default=None,
					help="""Path of databases. If dbPath isn't specified but a path is provided in "--database" option, this path of database will also be used in dbPath. 
					Otherwise, the program will search "database/" in program directory.
					[default: database/]""")

	p.add_argument( '-po','--pathogenOnly', action="store_true", required=False,
					help="Display pathogen only [default: None]")

	p.add_argument( '-da','--displayAll', action="store_true", required=False,
					help="Display all taxonomies including being filtered out [default: None]")

	args_parsed = p.parse_args()

	if not args_parsed.dbPath:
		if args_parsed.database and "/" in args_parsed.database:
			db_dir = search( '^(.*?)[^\/]+$', args_parsed.database )
			args_parsed.dbPath = db_dir.group(1)
		else:
			bin_dir = os.path.dirname(os.path.realpath(__file__))
			args_parsed.dbPath = bin_dir + "/database"

	return args_parsed

if __name__ == '__main__':
	argvs = parse_params()
	t.loadTaxonomy( argvs.dbPath )
	
	for line in sys.stdin:
		# skip header
		if line.startswith("LEVEL"):
			continue
		temp = line.split('\t')
		# only process species level
		if temp[0] != "species":
			continue
		# only process unfiltered
		if not argvs.displayAll and temp[15]:
			continue
		# display pathogen only
		if not temp[10] and argvs.pathogenOnly:
			continue

		lineage = t.taxid2lineage(temp[2])

		print( "%s\t%s" %
			( temp[13],
			  '\t'.join( lineage.split('|') )
		))