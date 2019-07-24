#!/usr/bin/env python

# Po-E (Paul) Li
# B-10, Los Alamos National Lab
# Date: 12/01/2016

import os
import sys

libPath = os.path.dirname(os.path.realpath(__file__))
taxonomyDir = libPath + "/database"
DEBUG=0

def loadPathogen(pathogen_file=taxonomyDir+"/pathogen.tsv"):
	if DEBUG: sys.stderr.write( "[INFO] Open pathogen file: %s\n" % pathogen_file )

	pathMeta = {}

	try:
		with open(pathogen_file, 'r', encoding='UTF-8') as f:
			for line in f:
				tid, name, source, location, host, disease = line.rstrip('\r\n').split('\t')
				pathMeta[tid]={}
				pathMeta[tid]['name'] = name
				pathMeta[tid]['source'] = source
				pathMeta[tid]['location'] = location
				pathMeta[tid]['host'] = host
				pathMeta[tid]['disease'] = disease
			f.close()
	except IOError:
		sys.stderr.write( "Failed to open pathogen file: %s.\n" % pathogen_file )

	if DEBUG: sys.stderr.write( "[INFO] Done parsing pathogen.tsv (%d pathogen loaded)\n" % len(pathMeta) )

	return pathMeta


if __name__ == '__main__':
	print(loadPathogen())