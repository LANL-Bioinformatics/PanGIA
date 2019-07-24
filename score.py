#!/usr/bin/env python3
import os
import sys
import pandas as pd
import taxonomy as t
from multiprocessing import Pool
__version__ = "1.0.0-RC5"

def loadUniqueness(uniqueness_file):
	uniq_info = {}

	try:
		with open(uniqueness_file, 'r', encoding='UTF-8') as f:
			for line in f:
				(tid, raw_len, sk, p, c, o, f, g, s, sn) = line.rstrip('\r\n').split('\t')
				
				if not tid in uniq_info: uniq_info[tid]={}
				
				uniq_info[tid]["raw"]          = int(raw_len)
				uniq_info[tid]["superkingdom"] = int(sk)
				uniq_info[tid]["phylum"]       = int(p) if p != '-' else p
				uniq_info[tid]["class"]        = int(c) if c != '-' else c
				uniq_info[tid]["order"]        = int(o) if o != '-' else o
				uniq_info[tid]["family"]       = int(f) if f != '-' else f
				uniq_info[tid]["genus"]        = int(g) if g != '-' else g
				uniq_info[tid]["species"]      = int(s) if s != '-' else s
				uniq_info[tid]["strain"]       = int(sn) if sn != '-' else sn
	except IOError:
		sys.stderr.write( "Failed to open uniqueness file: %s.\n" % uniqueness_file )

	return uniq_info

def scoreStrainUniqueness(res, tid, uniq_meta, method="default"):
	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']
	ranks_priority = dict(zip(ranks, range(8)))
	inlvl = res['LVL']
	p = 1
	u_len = 0
	t_len = 0

	if tid in uniq_meta:
		t_len = uniq_meta[tid]['raw']

	# add up score
	tol_score = 0

	#tid found in uniq metrics
	p_of_this_rank=1
	if tid in uniq_meta:
		(mapped, nm, read_len, read_num, read_rnr, dcnr) = (0,0,0,0,0,0)
		# interating all ranks from strain to superkingdom
		for lcr_lvl in ranks:
			# only procede when uniq info available in current rank 
			if uniq_meta[tid][lcr_lvl] != "-":
				u_len = uniq_meta[tid][lcr_lvl]
				u_pct = (t_len-u_len)/t_len
				
				#RC5
				p = 0.90*u_pct+0.05
				# v2.4.9
				# p = (u_pct+1)/2 #default

				if ":" in method:
					(low, hi) = method.split(':')
					try:
						low = float(low)
						hi = float(hi)

						if low > 1 and hi > 1:
							low = low/100
							hi = hi/100

						p = (hi-low)*u_pct+low
					except:
						sys.stderr.write( f"[ERROR] {method} is an invalid normalization parameter for scoring function.\n")
						sys.exit()
				else:
					if p>0.995: p=0.995
					if p<0.005: p=0.005

			# P(Gr) = Pr
			# P(Xr|Gr) = Cr
			# P(X) = 𝑃(𝑋_𝑟│𝐺_𝑟) ∙ 𝑃(𝐺_𝑟) + 𝑃(𝑋_𝑟│~𝐺_𝑟) ∙ 𝑃(~𝐺_𝑟)
			# P(G|X) = ∑ ∝ P(Xr|Gr)*P(G) / P(X)
			#        = ∑ ∝ ( C*P / (2CP-P-C+1) )

			# P(Gr) = Pr = p
			# prior propability of the genome presence in the sample

			if uniq_meta[tid][lcr_lvl] == "-": # the rank is more specific than target taxonomy
				if lcr_lvl in res:
					# add-up rank-specific info to target rank
					mapped   += res[lcr_lvl]["M"]
					nm       += res[lcr_lvl]["N"]
					read_len += res[lcr_lvl]["T"]
					read_num += res[lcr_lvl]["R"]
					read_rnr += res[lcr_lvl]["Rnr"]
					dcnr     += res[lcr_lvl]["DCnr"]
				else:
					continue
			else:
				if lcr_lvl in res or mapped>0:
					if lcr_lvl in res:
						mapped   += res[lcr_lvl]["M"]
						nm       += res[lcr_lvl]["N"]
						read_len += res[lcr_lvl]["T"]
						read_num += res[lcr_lvl]["R"]
						read_rnr += res[lcr_lvl]["Rnr"]
						dcnr     += res[lcr_lvl]["DCnr"]

					# skip this rank because the taxonomy is higher than this rank
					if uniq_meta[tid][lcr_lvl] == "-": continue

					# P(Xr|Gr) = Cr = c
					# Given genome exists, the probability of observing X number of reads mapped to this genome

					c = dcnr
					
					# if "propDCnr" in method:
					# 	propDCnr = dcnr/res["DCnr"]
					# 	c = propDCnr
					
					if c > 0.995:
						c = 0.995
					# else:
					# 	if c > 0.99:
					# 		c = 0.99

					try:
						#tol_score += (c*p/(2*c*p-p-c+1)) * read_rnr
						tol_score += (c*p/(2*c*p-p-c+1)) * p_of_this_rank
						p_of_this_rank *= 1-p
						(mapped, nm, read_len, read_num, read_rnr, dcnr) = (0,0,0,0,0,0)
					except:
						e = sys.exc_info()[0]
						sys.stderr.write( "[ERROR] %s.\n" % e )
						sys.stderr.write( "[INFO] t_len: %s, u_len: %s, lcr_lvl: %s.\n" % (t_len, u_len, lcr_lvl) )
						#sys.stderr.write( "[INFO] tid: %s, c: %s, p: %s, read_rnr: %s.\n" % (tid, c, p, read_rnr) )
						sys.exit()
				else:
					p_of_this_rank *= 1-p

		# P(G|X)
		#return tol_score/res["MR"] if tol_score/res["MR"] < 1 else 1
		if res["MRNr"] == 0:
			return 0
		else:
			#return tol_score/res["MRNr"] if tol_score/res["MRNr"] < 1 else 1
			return tol_score if tol_score < 1 else 1
	else:
		# Not in uniqueness matrix
		return "NA"

def clacRefMask(res_rollup, refs, bg_mask, tempdir, mb, mr, ml, mc, verbose, debug):
	depth={}
	# retrieve depths
	for ref in refs:
		# skip this ref if its organism doesn't pass thresholds
		(acc,slen,tid,tag) = ref.split('|')
		if not tid in res_rollup: continue
		if mb > res_rollup[tid]["rsMRNb"]: continue
		if mr > res_rollup[tid]["MR"]: continue
		if ml > res_rollup[tid]["LL"]: continue
		if mc > res_rollup[tid]["LC"]: continue

		depthfile = "%s/merged_sam/%s.sorted.depth"%(tempdir, ref)

		if not os.path.isfile(depthfile) or not os.path.getsize(depthfile):
			continue

		if ref in bg_mask:
			if debug: sys.stderr.write( "[DEBUG] Processing reference %s...\n"%ref)
			# load bg mask to dataframe
			maskbinstr = bin(bg_mask[ref])[2:]
			dfmask = pd.DataFrame( list(maskbinstr) )
			dfmask.columns = ['bg_mapped']
			dfmask.index += (int(slen)-len(maskbinstr)+1) 

			if debug: sys.stderr.write( "[DEBUG] Printing mask...\n")
			if debug: sys.stderr.write( str(dfmask)+"\n" )

			# load depth
			df = pd.read_csv(depthfile,
							sep='\t',
							names=['ref','pos','dep'],
							usecols=['pos','dep'],
							dtype={'pos':int,'dep':int},
							index_col=['pos']
			)

			if debug: sys.stderr.write( "[DEBUG] Printing base-by-base depths in target dataset...\n")
			if debug: sys.stderr.write( str(df)+"\n" )

			# join bg mask and the deps
			df1 = df.join(dfmask, on=df.index, how='left')
			df1['dep_mask'] = df1['dep']
			df1.loc[df1.bg_mapped=='1','dep_mask']=0

			if debug: sys.stderr.write( "[DEBUG] Printing MASKED base-by-base depths...\n")
			if debug: sys.stderr.write( str(df1)+"\n" )

			if tid in depth:
				depth[tid]['orig'].append(df1['dep'])
				depth[tid]['mask'].append(df1['dep_mask'])
			else:
				depth[tid]={}
				depth[tid]['orig'] = df1['dep'].copy()
				depth[tid]['mask'] = df1['dep_mask'].copy()

	return depth

def scoreBg(res_rollup, refs, bg_mask, tempdir, numthreads, mb, mr, ml, mc, verbose, debug, method="overlapping_prop_raw"):
	"""
	  Calculate a score by performing the 2-sample Kolmogorov-Smirnov test between 
	the depth of coverage (DoC) and background masked DoC.
	"""
	from scipy import stats
	import numpy as np
	if debug: sys.stderr.write( "[DEBUG] Using score.py version:%s\n"%__version__ )

	pool = Pool(processes=numthreads)
	jobs = []
	depth={}

	n = 10 #number of refs per chunk
	chunks = [refs[i:i+n] for i in range(0, len(refs), n)]

	for chunk in chunks:
		jobs.append( pool.apply_async(clacRefMask, (res_rollup, chunk, bg_mask, tempdir, mb, mr, ml, mc, verbose, debug)) )

	tol_jobs = len(jobs)
	cnt=0
	for job in jobs:
		depth_chunk = job.get()
		if depth_chunk:
			for tid in depth_chunk:

				lng = t.taxid2lineageDICT(tid, 1, 1)
				for rank in lng:
					if rank == "type": continue #skip "type" level
					#res_tree[pid][tid] = 1
					taxid = lng[rank]['taxid']

					# if no particular rank for this organism, use strain
					if taxid == 0:
						taxid = lng[rank]['name'] # for example: 'Flaviviridae - no_o_rank - no_c_rank'

					if taxid in depth:
						depth[taxid]['orig'].append(depth_chunk[tid]['orig'])
						depth[taxid]['mask'].append(depth_chunk[tid]['mask'])
					else:
						depth[taxid] = {}
						depth[taxid]['orig'] = depth_chunk[tid]['orig'].copy()
						depth[taxid]['mask'] = depth_chunk[tid]['mask'].copy()
		cnt+=1
		if verbose: print( "[INFO] Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100))

	#clean up
	pool.close()

	for tid in res_rollup:
		# skipping un-qualified taxas
		if mb > res_rollup[tid]["rsMRNb"] or mr > res_rollup[tid]["MR"] or ml > res_rollup[tid]["LL"] or mc > res_rollup[tid]["LC"]:
			continue

		if tid in depth:
			if debug: sys.stderr.write( "[DEBUG] Processing %s (taxid:%s)...\n"%(tid if '-' in tid else t.taxid2name(tid), tid))

			def transform_freq(orig, mask):
				dep = pd.DataFrame(data={
					'orig': orig, 'mask': mask
				})
				dep_co = dep.groupby(['orig']).count()
				dep_cm = dep.groupby(['mask']).count()
				dep_co.index.names = ['depth']
				dep_cm.index.names = ['depth']
				dep_comb = dep_co.join(dep_cm, how='outer').fillna(0).rename({'orig':'mask','mask':'orig'}, axis='columns')
				if debug: sys.stderr.write( "[DEBUG] Freq of binned depths...\n" )
				if debug: sys.stderr.write( str(dep_comb)+"\n" )
				return dep_comb
		
			def overlapping_prop_raw(orig, mask):
				non_overlapping_area = np.absolute(orig-mask).sum()
				sum_area = (orig+mask).sum()
				union_area = (sum_area+non_overlapping_area)/2
				non_overlapping_prop = non_overlapping_area/union_area
				if debug: sys.stderr.write( "[DEBUG] (non_overlapping_area, union_area, non_overlapping_prop) = (%s, %s, %s)\n"%(non_overlapping_area, union_area, non_overlapping_prop))
				return 1-non_overlapping_prop

			def overlapping_prop_freq(orig, mask):
				dep_comb = transform_freq(orig, mask)
				return overlapping_prop_raw(dep_comb['orig'], dep_comb['mask'])

			def KS_2samp_freq(orig, mask):
				"""
				2-sample K-S test
				"""
				dep_comb = transform_freq(orig, mask)
				d, pval = stats.ks_2samp( dep_comb['orig'], dep_comb['mask'] )
				if debug: sys.stderr.write( "[DEBUG] (d, pval)=(%s, %s)\n"%(d, pval) )
				return pval

			def chi2_contingency_freq(orig, mask):
				"""
				Chi2 test
				"""
				dep_comb = transform_freq(orig, mask)
				chi2, pval, dof, ex = stats.chi2_contingency([dep_comb['orig'], dep_comb['mask']])
				return pval

			def KS_2samp_freq_removedZero(orig, mask):
				dep_comb = transform_freq(orig, mask)
				if 0 in dep_comb.index:
					dep_comb = dep_comb.drop(0)
				d, pval = stats.ks_2samp( dep_comb['orig'], dep_comb['mask'] )
				if debug: sys.stderr.write( "[DEBUG] (d, pval)=(%s, %s)\n"%(d, pval) )
				return pval

			if debug: sys.stderr.write( "[DEBUG] Using method: %s...\n"%method)
			
			if method=="overlapping_prop_freq":
				res_rollup[tid]["S_BG"] = overlapping_prop_freq( depth[tid]['orig'], depth[tid]['mask'] )
			elif method=="KS_2samp_freq":
				res_rollup[tid]["S_BG"] = KS_2samp_freq( depth[tid]['orig'], depth[tid]['mask'] )
			elif method=="KS_2samp_freq_removedZero":
				res_rollup[tid]["S_BG"] = KS_2samp_freq_removedZero( depth[tid]['orig'], depth[tid]['mask'] )
			else:
				res_rollup[tid]["S_BG"] = overlapping_prop_raw( depth[tid]['orig'], depth[tid]['mask'] )

		else:
			res_rollup[tid]["S_BG"] = 1

	return res_rollup