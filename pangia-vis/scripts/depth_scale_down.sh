set -e

MERGE_SAM_DIR=$1
OUTDIR=$2

SCRIPTNAME=$0
SCRIPTNAME=$(dirname "$0")
BIN="$SCRIPTNAME/depth_scale_down.py"

mkdir -p $OUTDIR

# uniq organism tid
find $MERGE_SAM_DIR -type f -name '*.depth' | cut -d"|" -f3 | sort | uniq > $OUTDIR/unique_taxid.txt
parallel "$BIN $MERGE_SAM_DIR/*\|{}\|*.depth > $OUTDIR/{}.depth.scaledown" :::: $OUTDIR/unique_taxid.txt

set +e
