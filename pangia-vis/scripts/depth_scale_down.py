#!/usr/bin/env python3
import sys
import os
import glob
import subprocess

depth_files = sys.argv[1:]
ds_size = 2000
out_buffer = ""
block_size = 1

# read all lines from depth files
lines = []
for fn in depth_files:
    if os.path.isfile(fn):
        with open(fn) as f:
            lines += f.readlines()
        f.close()

# if total number of lines > 5000, downsale total number of lines to 5000
if len(lines) > ds_size:
    block_size = int(len(lines)/ds_size)

# parse lines
block_ref = None
block_start = 0
block_sum = 0
block_offset = 0
pos = 0
dep = 0
temp = []

if lines:
    for line in lines:
        if not line: next
        temp = line.strip().split('\t')
        if not len(temp)==3: next
        ref = temp[0]
        pos = int(temp[1])
        dep = int(temp[2])

        if block_ref != ref:
            if block_ref:
                block_offset += int(block_ref.split('|')[1])
            block_ref = ref
            block_start = 0

        if block_start==0:
            block_start = pos + block_offset
            block_sum = dep
        elif (block_start+block_size-1)<(pos+block_offset):
            out_buffer += "%s\t%s\t%d\n"%(temp[0], block_start, block_sum/block_size)
            block_start = pos + block_offset
            block_sum = dep
        else:
            block_sum += dep
    #last block
    out_buffer += "%s\t%s\t%d\n"%(temp[0], block_start, block_sum/(pos+block_offset-block_start+1))

# print to stdout
sys.stdout.write(out_buffer)
