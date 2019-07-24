var data = source.data
var filetext = 'NAME,TAXID,READ_COUNT,READ_COUNT_RNR,READ_COUNT_RSNB,PRI_READ_COUNT,REL_ABUNDANCE,SCORE,RPKM,LINEAR_LENGTH,LINEAR_COV,RS_DEPTH_COV_NR,DEPTH_COV,PATHOGEN,HOST,SOURCE,LOCATION,DISEASE\n'

for (i=0; i < data['taxa'].length; i++) {
    var currRow = [data['taxa'][i].toString(),
                    data['taxid'][i].toString(),
                    data['read'][i].toString(),
                    data['read_rnr'][i].toString(),
                    data['read_rsnb'][i].toString(),
                    data['read_pri'][i].toString(),
                    //data['reb_read'][i] === null ? 'NA' : data['reb_read'][i].toString(),
                    //data['reb_read_rnr'][i] === null ? 'NA' : data['reb_read_rnr'][i].toString(),
                    data['ra'][i].toString(),
                    data['score'][i].toString(),
                    data['rpkm'][i].toString(),
                    data['lnr_len'][i].toString(),
                    data['lnr_cov'][i].toString(),
                    data['rs_doc_nr'][i].toString(),
                    data['doc'][i].toString(),
                    data['pathogen'][i].toString(),
                    data['p_host'][i].toString(),
                    data['p_src'][i].toString(),
                    data['p_loc'][i].toString(),
                    data['p_dse'][i].toString().concat('\n')
    ]
                    
    var joined = currRow.join()
    filetext = filetext.concat(joined)
}

console.log(filetext)

var filename = 'pangia-vis_result.csv';
var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
}
else {
    var link = document.createElement("a")
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = "_blank"
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}