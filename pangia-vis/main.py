__name__      = "PanGIA-VIS: PanGIA result visualization tool"
__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = "1.0.0-RC3"
__date__      = "2018/09/10"
__copyright__ = "BSDv3"

import pandas as pd
import numpy as np
import sys
import os
from math import pi
from operator import itemgetter
import re
from os.path import dirname, join, isfile
from bokeh.plotting import curdoc, figure
from bokeh.layouts import row, layout, widgetbox, column
from bokeh.models import ColumnDataSource, HoverTool, Div, FactorRange, Range1d, TapTool, CustomJS
from bokeh.models.widgets import Button, Slider, Select, TextInput, RadioButtonGroup, DataTable, TableColumn, NumberFormatter, Panel, Tabs
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import RdYlBu11, Spectral4, Spectral11, Set1
from bokeh.util import logconfig

###############################################################################
# default cutoffs for dot-plot display and init variables
###############################################################################
def_val_patho      = 0     #display mode: 0 -> pathogen only, 1 -> all
def_val_min_len    = 50    #Minimum linear length [0-500] step=1
def_val_min_cov    = 0.01  #Minimum genome coverage [0-1] step=0.01
def_val_max_r_raw  = 100   #Minimum reads [0-500] step=1
def_val_max_r_rsnb = 10    #Minimum reads [0-100] step=1
def_val_min_score  = 0.5   #Minimum score [0-1] step=0.1
def_val_min_dc     = 10    #Minimum depth coverage MilliX (1X=1000mX) [0-10000] step=10
def_val_min_rsdc   = 1     #Minimum rank specific depth coverage in MilliX[0-1000] step=0.001

COLORS = RdYlBu11 #RdYlBu11
C_SET1 = Set1[9]
MIN_SIZE = 10
MAX_SIZE = 36

GCOV_CACHE={}
BSAT_CACHE={}
MASK_CACHE={}

df_filtered = pd.DataFrame()
bgfiles = []
bg_mask = {}

###############################################################################
# get pangia result
###############################################################################

args = curdoc().session_context.request.arguments
r_file = ""
proj_name = ""
r_basepath = ""
db_path = ""

if len(sys.argv) > 1:
    r_file = sys.argv[1]
    if isfile(r_file):
        logconfig.bokeh_logger.debug("[DEBUG] [INPUT] Input file by sys argv: %s"%r_file)
    else:
        exit("[ERROR] Input file by sys argv: not found.")

    try:
        r_file_fn = r_file.split('/')[-1].split('.')[0]
        proj_name = "(%s)"%r_file_fn
        r_basepath = "/".join( r_file.split("/")[:-1] )
    except:
        logconfig.bokeh_logger.debug("[DEBUG] [INPUT] Not able to retrieve file name: %s"%r_file)

if not r_file:
    try:
        r_file = args.get('r')[0].decode("utf-8")
        if isfile(r_file):
            logconfig.bokeh_logger.debug("[DEBUG] [INPUT] Input file by URL: %s"%r_file)
            r_basepath = "/".join( r_file.split("/")[:-1] )
        else:
            exit("[ERROR] Input file by URL: not found.")
        
        try:
            proj_name = "(%s)"%args.get('p')[0].decode("utf-8")
        except:
            logconfig.bokeh_logger.debug("[DEBUG] [INPUT] No project name found (p).")
    except:
        exit("[ERROR] No PanGIA result provided.")

# get file path + prefix
r_path = r_file.split('.')[0]

# load data and cleanup
res_df = pd.read_csv(r_file, sep='\t')
res_df = res_df.replace('NA', np.nan)
res_df = res_df.fillna(0)
res_df['NAME'] = res_df['NAME'].str.replace(":", " ")
res_df['PATHOGEN'] = res_df['PATHOGEN'].apply(str)

logconfig.bokeh_logger.debug("[DEBUG] [INPUT] Done loading input file.")

###############################################################################
# dot plot
###############################################################################

ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

cols_map = {
    'None': 'None',
    'Read count (raw)': 'READ_COUNT',
    'Read count (raw, rebalanced)': 'REB_READ_COUNT',
    'Read count (raw, only primaly align)': 'PRI_READ_COUNT',
    'Read count (RNR)': 'READ_COUNT_RNR',
    'Read count (RNR, rebalanced)': 'REB_READ_COUNT_RNR',
    'Read count (RSNB)': 'READ_COUNT_RSNB',
    'Genome coverage': 'LINEAR_COV',
    'RPKM': 'RPKM',
    'Depth of coverage (RSNR)': 'RS_DEPTH_COV_NR',
    'Depth of coverage': 'DEPTH_COV',
    'Score': 'SCORE',
    'Score (uniq)': 'SCORE_UNIQ',
    'Score (background)': 'SCORE_BG',
    'Relative abundance': 'REL_ABUNDANCE'
}

# plot chart
# plot DataSource
dp_source = ColumnDataSource(data=dict(
    taxa=[], y=[], color=[], size=[], read=[], read_rnr=[], read_pri=[], read_rsnb=[], score=[],
    rpkm=[], lnr_len=[], lnr_cov=[], rs_doc_nr=[], doc=[], ra=[],
    g_size=[], pathogen=[], disease=[], host=[], location=[],
    STR=[],SPE=[],GEN=[],FAM=[],ORD=[],CLA=[],PHY=[],SK=[],ROOT=[],
    STR_rnb=[],SPE_rnb=[],GEN_rnb=[],FAM_rnb=[],ORD_rnb=[],CLA_rnb=[],PHY_rnb=[],SK_rnb=[],ROOT_rnb=[],
    STR_rnr=[],SPE_rnr=[],GEN_rnr=[],FAM_rnr=[],ORD_rnr=[],CLA_rnr=[],PHY_rnr=[],SK_rnr=[],ROOT_rnr=[],
    gcov_ref=[], gcov_pos=[], gcov_dep=[], gcov_col=[],
    bsat_ref=[], bsat_str=[], bsat_end=[], bsat_col=[], best_not=[]
))

hover = HoverTool(tooltips=[
    ("Name", "@taxa"),
    ("TaxID", "@taxid"),
    ("Score", "@score{0.00}"),
    ("Linear Cov", "@lnr_cov{0,0.00}"),
    ("Depth", "@doc{0,0.00}"),
    ("Depth (RSNR)", "@rs_doc_nr{0,0.00}"),
    ("Relative abundance", "@ra{0,0.00}"),
    ("Read (raw)", "@read{0,0}"),
    ("Read (rnr)", "@read_rnr{0,0.00}"),
    ("Read (rsnb)", "@read_rsnb{0,0.00}"),
    ("Read (primary)", "@read_pri{0,0}"),
    ("RPKM", "@rpkm{0,0.00}"),
    ("Score (uniq)", "@score_uniq{0.00}"),
    ("Score (background)", "@score_bg{0.00}"),
    ("Genome Size (bp)", "@g_size{0,0}"),
    ("Pathogen", "@pathogen"),
])

p = figure(plot_height=800, plot_width=800, title="", y_axis_type="log",
    toolbar_location="below", toolbar_sticky=False, 
    x_range=FactorRange(), y_range=Range1d(),
    output_backend="webgl",
    tools=["pan,wheel_zoom,box_zoom,reset,tap",hover]
)

p.circle(x="taxa", y="y",
    source=dp_source, color="color",
    size="size", line_color='#C9C9C9', alpha=0.6, hover_alpha=0.5
)

p.toolbar.logo = None

###############################################################################
# Dashboard
###############################################################################

###### input reads distribution
pieHover = HoverTool(
    tooltips=[
        ( 'name', '@name' ),
        ( 'reads', '@val{,} (@pct{%0.0f})' ),
    ])

pieInReadsDS = ColumnDataSource(data=dict( name=['NA'], start_angle=[0], end_angle=[2*pi], color=['#EFF0F1'], val=['NA'], pct=['NA']) )
pieInReadsFigure = figure( x_range=(-1.3, 4), output_backend="webgl", y_range=(-2, 2), plot_width=380, plot_height=220, title="Total reads:", tools=[pieHover] )

pieInReadsFigure.annular_wedge(
    x=0, y=0, alpha=0.7,
    legend='name', start_angle='start_angle', end_angle='end_angle', color='color',
    inner_radius=0.7, outer_radius=1.2, source=pieInReadsDS)

###### Flag distribution
pieFlagDS = ColumnDataSource(data=dict( name=['NA'], start_angle=[0], end_angle=[2*pi], color=['#EFF0F1'], val=['NA'], pct=['NA']) )
pieFlagFigure = figure( x_range=(-1.5, 4), y_range=(-2, 2), plot_width=380, plot_height=220, title="Target reads distribution:", tools=[pieHover] )

pieFlagFigure.annular_wedge(
    x=0, y=0, alpha=0.7,
    legend='name', start_angle='start_angle', end_angle='end_angle', color='color',
    inner_radius=0.7, outer_radius=1.2, source=pieFlagDS)

###### pathogen stats
piePathoDS = ColumnDataSource(data=dict( name=['NA'], start_angle=[0], end_angle=[2*pi], color=['#EFF0F1'], val=['NA'], pct=['NA']) )
piePathoFigure = figure( x_range=(-1.27, 4), output_backend="webgl", y_range=(-2, 2), plot_width=380, plot_height=220, title="Pathogen reads distribution:", tools=[pieHover] )

piePathoFigure.annular_wedge(
    x=0, y=0, alpha=0.7,
    legend='name', start_angle='start_angle', end_angle='end_angle', color='color',
    inner_radius=0.7, outer_radius=1.2, source=piePathoDS
)

###############################################################################
# data table
###############################################################################

table_cols = [
    TableColumn(field="taxa",     title="Name", width=800),
    TableColumn(field="taxid",    title="TaxID"),
    TableColumn(field="read_rnr", title="Read (RNR)",    formatter=NumberFormatter(format='0,0.00')),
    TableColumn(field="read_rsnb",title="Read (RSNB)",   formatter=NumberFormatter(format='0,0.00')),
    TableColumn(field="read_pri", title="Read (pri)",    formatter=NumberFormatter(format='0,0')),
    TableColumn(field="score",    title="Score",         formatter=NumberFormatter(format='0.00')),
    TableColumn(field="score_uniq",title="Score (uniq)", formatter=NumberFormatter(format='0.00')),
    TableColumn(field="score_bg", title="Score (bg)",    formatter=NumberFormatter(format='0.00')),
    TableColumn(field="rpkm",     title="RPKM",          formatter=NumberFormatter(format='0,0')),
    TableColumn(field="lnr_cov",  title="Genome Cov",    formatter=NumberFormatter(format='0,0.00')),
    TableColumn(field="rs_doc_nr",title="Depth (RSNR)",  formatter=NumberFormatter(format='0,0.00')),
    TableColumn(field="doc",      title="Depth",         formatter=NumberFormatter(format='0,0.00')),
    TableColumn(field="ra",       title="Rel Abundance", formatter=NumberFormatter(format='0.00%')),
    TableColumn(field="pathogen", title="Pathogen"),
]

#use "index_position=None" instead of "row_headers=False" on bokeh ver >0.12.10
data_table = DataTable(source=dp_source, columns=table_cols, index_position=None, width=1070, height=200)
result_table = widgetbox(data_table, width=1070)

###############################################################################
# Rank specific bar chat
###############################################################################

# Rank specific DataSource
rs_source = ColumnDataSource(data=dict(
    rs_rank = ["SKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"],
    rs_raw =[0,0,0,0,0,0,0,0],
    rs_rnr =[0,0,0,0,0,0,0,0],
    rs_rnb =[0,0,0,0,0,0,0,0],
    rs_ri  =[0,0,0,0,0,0,0,0],
))

hover_read = HoverTool(tooltips=[
    ("Read(raw)", "@rs_raw{0,0}"),
    ("Read(rnr)", "@rs_rnr{0,0.00}"),
    ("Read(rsnb)", "@rs_rnb{0,0.00}"),
])

hover_read_idt = HoverTool(tooltips=[
    ("identity", "@rs_ri{0,0.00%}"),
])

bar_opts = dict(line_color="white", height=0.48)

rs_p = figure(plot_width=270, plot_height=210,
    x_axis_type="log",
    title="Rank specific read counts",
    y_range=rs_source.data['rs_rank'],
    x_range=(0.1,600), 
    toolbar_location=None,
    output_backend="webgl",
    tools=[hover_read]
)

rs_p.hbar(y="rs_rank", left=0.1,      right='rs_rnb', color=COLORS[9], source=rs_source, **bar_opts)
rs_p.hbar(y="rs_rank", left='rs_rnb', right='rs_rnr', color=COLORS[6], source=rs_source, **bar_opts)
rs_p.hbar(y="rs_rank", left='rs_rnr', right='rs_raw', color=COLORS[2], source=rs_source, **bar_opts)
rs_p.xaxis.major_label_orientation = 3.1415926/4
rs_p.toolbar.logo = None

rsi_p = figure(plot_width=270, plot_height=210,
    #x_axis_type="log",
    title="Average mapping identity",
    y_range=rs_source.data['rs_rank'],
    x_range=(0,1),
    toolbar_location=None,
    output_backend="webgl",
    tools=[hover_read_idt]
)
rsi_p.hbar(y="rs_rank", left=0, right='rs_ri', color=COLORS[9], source=rs_source, **bar_opts)
rsi_p.xaxis.major_label_orientation = 3.1415926/4

###############################################################################
# Genome Depth & BSAT
###############################################################################

# Genome coverage DataSource
hover_gcov = HoverTool(tooltips=[
    ("Depth","information"),
    ("Reference", "@gc_ref"),
    ("Position",  "@gc_pos"),
    ("Depth",     "@gc_dep")],
    mode='vline',
    names=['b_gcov']
)

hover_bsat = HoverTool(tooltips=[
    ("Marker","information"),
    ("Reference", "@b_ref"),
    ("Location",  "@b_str..@b_end"),
    ("Note",      "@b_not")],
    mode='vline',
    names=['b_bsat']
)

hover_mask = HoverTool(tooltips=[
    ("Background","mask information"),
    ("Reference", "@m_ref"),
    ("Location",  "@m_str..@m_end")],
    names=['b_mask']
)

gcov_source = ColumnDataSource(data=dict(
    gc_ref = [],
    gc_pos = [],
    gc_dep = [],
    gc_col = []
))

bsat_source = ColumnDataSource(data=dict(
    b_ref = [],
    b_str = [],
    b_end = [],
    b_col = [],
    b_not = [],
))

mask_source = ColumnDataSource(data=dict(
    m_ref = [],
    m_str = [],
    m_end = []
))

coverage_p_log = figure(
    plot_width=1060,
    plot_height=190,
    min_border=0,
    y_axis_label='Depth (x)',
    #x_range=(1, 1000),
    #y_range=(0.1, 100),
    output_backend="webgl",
    y_axis_type="log",
    title='Not available at this rank (strain only)',
    tools=["wheel_zoom,box_zoom,pan,reset,save",hover_gcov,hover_bsat,hover_mask]
)

coverage_p_lin = figure(
    plot_width=1060,
    plot_height=190,
    min_border=0,
    y_axis_label='Depth (x)',
    #x_range=(1, 1000),
    #y_range=(0.1, 100),
    output_backend="webgl",
    title='Not available at this rank (strain only)',
    tools=["wheel_zoom,box_zoom,pan,reset,save",hover_gcov,hover_bsat,hover_mask]
)

# coverage
coverage_p_log.vbar(x='gc_pos', top='gc_dep', width=1, color='gc_col', bottom=0.1, alpha=0.2, name="b_gcov", source=gcov_source)
coverage_p_lin.vbar(x='gc_pos', top='gc_dep', width=1, color='gc_col', bottom=0,   alpha=0.2, name="b_gcov", source=gcov_source)
coverage_p_log.y_range.start = 0.1
coverage_p_lin.y_range.start = 0
coverage_p_log.x_range.start = 0
coverage_p_lin.x_range.start = 0

# additional y-range
coverage_p_log.extra_y_ranges = {"BSAT": Range1d(start=0.1, end=10), "MASK": Range1d(start=0.1, end=10)}
coverage_p_lin.extra_y_ranges = {"BSAT": Range1d(start=0,   end=10), "MASK": Range1d(start=0, end=10)}

# BSAT
coverage_p_log.hbar( left='b_str', right='b_end', color='b_col', y=0.3, height=0.13, source=bsat_source, name="b_bsat", y_range_name="BSAT")
coverage_p_lin.hbar( left='b_str', right='b_end', color='b_col', y=2.3, height=1.0,  source=bsat_source, name="b_bsat", y_range_name="BSAT")

# background mask
coverage_p_log.quad( top=10, bottom=0.1, left='m_str', right='m_end', color='#D3D3D3', line_width=0, fill_alpha=0.8, source=mask_source, name="b_mask", y_range_name="MASK")
coverage_p_lin.quad( top=10, bottom=0.1, left='m_str', right='m_end', color='#D3D3D3', line_width=0, fill_alpha=0.8, source=mask_source, name="b_mask", y_range_name="MASK")

coverage_p_log.toolbar.logo = None
coverage_p_lin.toolbar.logo = None
tab1 = Panel(child=coverage_p_log, title="Log")
tab2 = Panel(child=coverage_p_lin, title="Linear")
coverage_p_tabs = Tabs(tabs=[tab1, tab2])

###############################################################################
# Functions
###############################################################################

def parseLog( logfile ):
    total_input=0
    total_mapped=0
    info={}
    info['Target']=0
    info['Host']=0
    info['Ignored']=0

    with open(logfile) as f:
        for line in f:
            if "Total number of input reads" in line:
                reg = re.search( 'Total number of input reads: (\d+)', line )
                total_input = int(reg.group(1))
            elif "Total number of mapped reads" in line:
                reg = re.search( 'Total number of mapped reads: (\d+)', line )
                total_mapped = int(reg.group(1))
            elif "Total number of host reads" in line:
                reg = re.search( 'Total number of host reads: (\d+)', line )
                info['Host'] = int(reg.group(1))
            elif "Done processing SAM file" in line:
                reg = re.search( 'Done processing SAM file, (\d+) alignment', line )
                total_mapped = int(reg.group(1))
            elif "Input background  : [" in line:
                reg = re.search( 'Input background  : \[(.+)\]', line )
                bg_file_string = reg.group(1)
                bg_file_string = bg_file_string.strip("'")
                global bgfiles
                bgfiles = bg_file_string.split("', '")
            elif "Database          : [" in line:
                reg = re.search( 'Database          : \[(.+)\]', line )
                db_string = reg.group(1)
                db_string = db_string.strip("'")
                global db_path
                db_file = db_string.split("', '")[0]
                db_path = "/".join( db_file.split("/")[:-1] )
            elif "Total number of ignored reads" in line:
                reg = re.search( 'Total number of ignored reads .*: (\d+)', line )
                info['Ignored'] = int(reg.group(1))
                break
    
    info['Target'] = total_mapped - info['Host'] - info['Ignored']

    if total_input > 0:
        info['Unmapped'] = total_input - total_mapped

    return info

def genPieValues( distinfo ):
    tol_cnt = 0
    ang_start=0
    ang_size=6.27
    colors= list(itemgetter(9, 3, 0, 1, 2, 5, 7, 8, 9)(Spectral11))
    
    p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct = [],[],[],[],[],[]
    
    for name in distinfo:
        p_val.append(distinfo[name])
        tol_cnt += distinfo[name]
    
    percents=[0]
    for name in distinfo:
        pct=distinfo[name]/tol_cnt
        p_pct.append(pct)
        p_name.append( "%.1f%% %s"%(pct*100,name))
        percents.append( percents[-1]+(distinfo[name]/tol_cnt) )
        
    p_str_ang = [p*2*pi for p in percents[:-1]]
    p_stp_ang = [p*2*pi for p in percents[1:]]
    
    p_col = colors[:len(p_name)]
    
    return (p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct)

def setupPieFigures( pieFigure ):
    pieFigure.axis.visible = False
    pieFigure.grid.visible = False
    pieFigure.legend.location = "center_right"
    pieFigure.toolbar.logo = None
    pieFigure.toolbar_location = None
    pieFigure.outline_line_width = 0
    pieFigure.outline_line_alpha = 0

def result_filter():
    selected = res_df
    if 'RS_DEPTH_COV_NR' in res_df: #backward compatibility
        selected = res_df[
            (res_df.LEVEL           == rank.value) &
            (res_df.READ_COUNT      >= max_r_raw.value) &
            (res_df.READ_COUNT_RSNB >= max_r_rsnb.value) &
            (res_df.LINEAR_COV      >= min_cov.value) &
            (res_df.SCORE           >= min_score.value) &
            (res_df.DEPTH_COV       >= min_dc.value/1000) &
            (res_df.RS_DEPTH_COV_NR >= min_rsdc.value/1000)
        ]
    else:
        selected = res_df[
            (res_df.LEVEL           == rank.value) &
            (res_df.READ_COUNT      >= max_r_raw.value) &
            (res_df.READ_COUNT_RSNB >= max_r_rsnb.value) &
            (res_df.LINEAR_COV      >= min_cov.value) &
            (res_df.SCORE           >= min_score.value) &
            (res_df.DEPTH_COV       >= min_dc.value/1000)
        ]

    # additional cutoff filed
    if add_cut_c.value != "None":
        col = cols_map[add_cut_c.value]
        if col in selected:
            selected = selected[ selected[col] > float(add_cut_v.value) ]
        else:
            logconfig.bokeh_logger.debug("[DEBUG] [FILTER] Column %s not found."%col)

    #if patho.active == 'Pathogen Only':
    if patho.active == 0:
        selected = selected[selected.PATHOGEN=='Pathogen']
    if (org_n.value != ""):
        selected = selected[selected.NAME.str.contains(org_n.value)==True]
    if (dse_n.value != ""):
        selected = selected[selected.DISEASE.str.contains(dse_n.value)==True]

    return selected

def update():
    df = result_filter()

    y_col = cols_map[y.value]
    c_col = cols_map[color.value]
    s_col = cols_map[size.value]

    p.x_range.factors = []
    p.title.text = "Loading..."

    if len(df):
        #init
        p.x_range.factors = sorted(df['NAME'].tolist())
        p.y_range.start = 1
        p.y_range.end = 2
        if len(df[y_col].values) > 0:
            p.y_range.end = max(df[y_col].values)*5   

        p.xaxis.axis_label = rank.value
        p.yaxis.axis_label = y.value
        p.xaxis.major_label_orientation = 3.1415926/4

        #colors
        c=[COLORS[8]]*len(df['NAME'].values)
        if color.value != 'None':
            c = [ COLORS[int(x*110/11)] for x in df[c_col].values ]

        #sizes
        sz=[10]*len(df['NAME'].values)
        if size.value != 'None':
            max_x = max(df[s_col].values) if max(df[s_col].values) > 0 else 0.001
            sz = [ int(MIN_SIZE + x/max_x*(MAX_SIZE-MIN_SIZE)) for x in df[s_col].values ]

        # patho stats for pie chart
        info = {}
        info['Pathogen'] = df.loc[df.PATHOGEN=="Pathogen","PRI_READ_COUNT"].sum()
        info['Not pathogen'] = df.loc[:,"PRI_READ_COUNT"].sum() - info['Pathogen']
        (p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct) = genPieValues( info )
        piePathoDS.data=dict(
            name=p_name, start_angle=p_str_ang, end_angle=p_stp_ang, color=p_col, val=p_val, pct=p_pct
        )

        # target stats for pie chart
        if "FLAG" in df:
            info = {}
            for flag in df.FLAG.unique():
                if flag == 'B':
                    name = "Bacteria"
                elif flag == 'A':
                    name = "Archaea"
                elif flag == 'E':
                    name = "Eukaryota"
                elif flag == 'V':
                    name = "Viruses"
                else:
                    name = flag
                
                info[name] = df.loc[df.FLAG==flag, "PRI_READ_COUNT"].sum()

            (p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct) = genPieValues( info )
            pieFlagDS.data=dict(
                name=p_name, start_angle=p_str_ang, end_angle=p_stp_ang, color=p_col, val=p_val, pct=p_pct
            )

        #set dp_source data
        p.title.text = "Rendering..."
        dp_source.data = dict(
            taxa      = df['NAME'].tolist(),
            y         = df[y_col].tolist(),
            color     = c,
            size      = sz,
            rank      = df['LEVEL'].tolist(),
            read      = df['READ_COUNT'].tolist(),
            read_rnr  = df['READ_COUNT_RNR'].tolist(),
            read_rsnb = df['READ_COUNT_RSNB'].tolist(),
            read_pri  = df['PRI_READ_COUNT'].tolist(),
            score     = df['SCORE'].tolist(),
            score_bg  = df['SCORE_BG'].tolist() if 'SCORE_BG' in df else [None]*len(c),
            score_uniq= df['SCORE_UNIQ'].tolist(),
            rpkm      = df['RPKM'].tolist(),
            lnr_len   = df['LINEAR_LENGTH'].tolist(),
            lnr_cov   = df['LINEAR_COV'].tolist(),
            rs_doc_nr = df['RS_DEPTH_COV_NR'].tolist() if 'RS_DEPTH_COV_NR' in df else [None]*len(c),
            doc       = df['DEPTH_COV'].tolist(),
            ra        = df['REL_ABUNDANCE'].tolist(),
            pathogen  = [ "Yes" if x=="Pathogen" else "No" for x in df['PATHOGEN'].values ],
            p_host    = df['HOST'].tolist(),
            p_src     = df['SOURCE'].tolist(),
            p_loc     = df['LOCATION'].tolist(),
            p_dse     = df['DISEASE'].tolist(),
            taxid     = [ str(x).replace(".0","") for x in df['TAXID'].values ],
            g_size    = df['TOL_GENOME_SIZE'].tolist(),
            STR       = df['STR'].tolist(),
            SPE       = df['SPE'].tolist(),
            GEN       = df['GEN'].tolist(),
            FAM       = df['FAM'].tolist(),
            ORD       = df['ORD'].tolist(),
            CLA       = df['CLA'].tolist(),
            PHY       = df['PHY'].tolist(),
            SK        = df['SK'].tolist(),
            ROOT      = df['ROOT'].tolist(),
            STR_rnb   = df['STR_rnb'].tolist(),
            SPE_rnb   = df['SPE_rnb'].tolist(),
            GEN_rnb   = df['GEN_rnb'].tolist(),
            FAM_rnb   = df['FAM_rnb'].tolist(),
            ORD_rnb   = df['ORD_rnb'].tolist(),
            CLA_rnb   = df['CLA_rnb'].tolist(),
            PHY_rnb   = df['PHY_rnb'].tolist(),
            SK_rnb    = df['SK_rnb'].tolist(),
            ROOT_rnb  = df['ROOT_rnb'].tolist(),
            STR_rnr   = df['STR_rnr'].tolist(),
            SPE_rnr   = df['SPE_rnr'].tolist(),
            GEN_rnr   = df['GEN_rnr'].tolist(),
            FAM_rnr   = df['FAM_rnr'].tolist(),
            ORD_rnr   = df['ORD_rnr'].tolist(),
            CLA_rnr   = df['CLA_rnr'].tolist(),
            PHY_rnr   = df['PHY_rnr'].tolist(),
            SK_rnr    = df['SK_rnr'].tolist(),
            ROOT_rnr  = df['ROOT_rnr'].tolist(),
            STR_ri   = df['STR_ri'].tolist() if 'STR_ri' in df else [None]*len(c),
            SPE_ri   = df['SPE_ri'].tolist() if 'SPE_ri' in df else [None]*len(c),
            GEN_ri   = df['GEN_ri'].tolist() if 'GEN_ri' in df else [None]*len(c),
            FAM_ri   = df['FAM_ri'].tolist() if 'FAM_ri' in df else [None]*len(c),
            ORD_ri   = df['ORD_ri'].tolist() if 'ORD_ri' in df else [None]*len(c),
            CLA_ri   = df['CLA_ri'].tolist() if 'CLA_ri' in df else [None]*len(c),
            PHY_ri   = df['PHY_ri'].tolist() if 'PHY_ri' in df else [None]*len(c),
            SK_ri    = df['SK_ri'].tolist()  if 'SK_ri' in df else [None]*len(c),
            ROOT_ri  = df['ROOT_ri'].tolist() if 'ROOT_ri' in df else [None]*len(c),
        )
    else:
        dp_source.data = dict(taxa=[],y=[],color=[],size=[])
    
    p.title.text = "%d taxonomies profiled"%len(df)

    # save filtered dataframe
    global df_filtered
    df_filtered = df

def cutoff_toggle():
    if max_r_raw.value == 0 and max_r_rsnb.value == 0 and min_cov.value == 0 and min_score.value == 0 and min_dc.value == 0 and min_rsdc.value == 0:
        min_cov.value = def_val_min_cov
        max_r_raw.value = def_val_max_r_raw
        max_r_rsnb.value = def_val_max_r_rsnb
        min_score.value = def_val_min_score
        min_dc.value = def_val_min_dc
        min_rsdc.value = def_val_min_rsdc
    else:
        max_r_raw.value = 0
        max_r_rsnb.value = 0
        min_cov.value = 0
        min_score.value = 0
        min_dc.value = 0
        min_rsdc.value = 0
    update()

def loading_log():
    # try loading log
    logfile = r_path + ".pangia.log"
    # skip empty files
    if os.path.isfile(logfile) and os.path.getsize(logfile) > 0:
        logconfig.bokeh_logger.debug("[DEBUG] [LOADING_LOG] Found logfile: %s"%logfile)

        info = parseLog( logfile )
        (p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct) = genPieValues( info )
        pieInReadsDS.data=dict(
            name=p_name,
            start_angle=p_str_ang,
            end_angle=p_stp_ang,
            color=p_col,
            val=p_val,
            pct=p_pct
        )

        tol_input_reads=0
        for x in info:
            tol_input_reads += info[x]
    else:
        pieInReadsFigure.title.text = "Total %s reads:"%tol_input_reads
        logconfig.bokeh_logger.debug("[DEBUG] [LOADING_LOG] The logfile not found.")

def loadBgMask(bgfiles):
    """
    Loading background mask files
    """
    import gzip
    import json

    logconfig.bokeh_logger.debug( "[DEBUG] [LOAD_BG] Input bgfiles: %s."%str(bgfiles) )

    global bg_mask
    for bgfile in bgfiles:
        if not os.path.isfile(bgfile):
            if os.path.isfile(r_basepath+"/"+bgfile):
                bgfile = r_basepath+"/"+bgfile
            else:
                continue
            logconfig.bokeh_logger.debug( "[DEBUG] [LOAD_BG] try loading: %s."%bgfile )

        try:
            with gzip.GzipFile(bgfile, 'r') as f:
                mask = json.loads(f.read())
                for ref in mask:
                    maskint = int(mask[ref][2:], 16)
                    if ref in bg_mask:
                        bg_mask[ref] |= maskint
                    else:
                        bg_mask[ref] = maskint
                logconfig.bokeh_logger.debug( "[DEBUG] [LOAD_BG] JSON file %s loaded."%bgfile )
        except IOError:
            logconfig.bokeh_logger.debug( "[DEBUG] [LOAD_BG] Failed to open background mask files: %s."%bgfile )

    logconfig.bokeh_logger.debug( "[DEBUG] [LOAD_BG] Done loading background JSON files." )

def update_rankspec_rc(attr, old, new):
    idx=0
    try:
        idx = new[0]
    except:
        pass

    rs_fields = {
        'rs_raw': ['SK','PHY','CLA','ORD','FAM','GEN','SPE','STR'],
        'rs_rnr': ['SK_rnr','PHY_rnr','CLA_rnr','ORD_rnr','FAM_rnr','GEN_rnr','SPE_rnr','STR_rnr'],
        'rs_rnb': ['SK_rnb','PHY_rnb','CLA_rnb','ORD_rnb','FAM_rnb','GEN_rnb','SPE_rnb','STR_rnb'],
        'rs_ri' : ['SK_ri','PHY_ri','CLA_ri','ORD_ri','FAM_ri','GEN_ri','SPE_ri','STR_ri']
    }

    rs_fields_vals = {
        'rs_raw': [],
        'rs_rnr': [],
        'rs_rnb': [],
        'rs_ri' : []
    }

    if len(dp_source.data['taxa'])>0:
        rs_p.background_fill_color = "white"
        rsi_p.background_fill_color = "white"
        for rs_type in rs_fields:
            for field in rs_fields[rs_type]:
                try:
                    rs_fields_vals[rs_type].append( dp_source.data[field][idx] )
                except:
                    rs_fields_vals[rs_type].append( 0 )

            rs_source.data[rs_type] = rs_fields_vals[rs_type]

        max_raw_cnt = 100
        try:
            max_raw_cnt = np.max(rs_source.data['rs_raw'])
        except:
            pass
        rs_p.x_range.end = max_raw_cnt*1.5

        try:
            ids = [x for x in rs_source.data['rs_ri'] if x > 0]
            min_id = np.min(ids)
            rsi_p.x_range.start = min_id-0.05
        except:
            rsi_p.background_fill_color = "#F3F3F3"
            pass
    else:
        rs_p.background_fill_color = "#F3F3F3"
        rsi_p.background_fill_color = "#F3F3F3"
        for rs_type in rs_fields:
            rs_source.data[rs_type] = [0]*len(rs_source.data[rs_type])

def update_genome_cov(attr, old, new):
    def reset_genome_plot(code=0):
        """
        Reset the genome plot and provide additional info
        """
        gcov_source.data = dict( gc_ref = [], gc_pos = [], gc_dep = [], gc_col = [])
        bsat_source.data = dict( b_ref = [], b_str = [], b_end = [], b_col = [], b_not = [])
        mask_source.data = dict( m_ref = [], m_str = [], m_end = [])
        text = "Not available"
        if code == 1:
            text = "Not available at this rank (strain only)"
        if code == 2:
            text = "Not available (unknown)"
        coverage_p_log.title.text = text
        coverage_p_lin.title.text = text
        coverage_p_log.background_fill_color = "#F3F3F3"
        coverage_p_lin.background_fill_color = "#F3F3F3"
        return

    if len(dp_source.data['taxa']):
        idx=0
        if len(new):
            try:
                idx = new[0]
            except:
                pass

        try:
            if dp_source.data['rank'][idx] != 'strain':
                reset_genome_plot(1)
                return
        except:
            reset_genome_plot()
            return

        # Genome coverage files
        name = dp_source.data['taxa'][idx]
        taxid = dp_source.data['taxid'][idx]
        gsize = dp_source.data['g_size'][idx]
        fn = r_path + "/" + taxid + ".depth.scaledown"
        pu = pd.DataFrame()
        mask_bp = 0

        # seting the end of x-axis
        logconfig.bokeh_logger.debug( f"[DEBUG] [GENO_COV] x_range.end: {gsize}" )
        coverage_p_log.x_range.end = gsize
        coverage_p_lin.x_range.end = gsize

        logconfig.bokeh_logger.debug( "[DEBUG] [GENO_COV] try loading depth file: %s."%fn )
        if os.path.isfile(fn) and os.path.getsize(fn) > 0:
            coverage_p_log.background_fill_color = "white"
            coverage_p_lin.background_fill_color = "white"
            # check if in the cache
            if taxid in GCOV_CACHE:
                gcov_source.data = GCOV_CACHE[taxid]
            else:
                coverage_p_log.title.text = f"Loading genome coverage for {name}..."
                coverage_p_lin.title.text = f"Loading genome coverage for {name}..."
                #loading coverage files
                pu = pd.read_csv( fn, sep='\t', header=None, names=['ref', 'pos', 'dep'] )
                pu = pu.drop_duplicates(subset=['ref', 'pos'])

                #ref colors
                ref_color = {}
                ref_names = pu['ref'].unique().tolist()
                COLORS14 = Spectral4[:3]+COLORS
                for idx,ref in enumerate(ref_names):
                    ref_color[ref] = COLORS14[idx%14]

                gcov_source.data = dict(
                    gc_ref = pu['ref'].tolist(),
                    gc_pos = pu['pos'].tolist(),
                    gc_dep = pu['dep'].tolist(),
                    gc_col = [ref_color[x] for x in pu['ref'].tolist()]
                )

                GCOV_CACHE[taxid] = gcov_source.data

            try:
                max_dep = max(gcov_source.data['gc_dep'])
                coverage_p_log.y_range.end = max_dep*3
                coverage_p_lin.y_range.end = max_dep*1.1
                logconfig.bokeh_logger.debug( f"[DEBUG] [GENO_COV] Max depth: {max_dep}" )
                logconfig.bokeh_logger.debug( f"[DEBUG] [GENO_COV] coverage_p_log.y_range.end: {coverage_p_log.y_range.end}" )
                logconfig.bokeh_logger.debug( f"[DEBUG] [GENO_COV] coverage_p_lin.y_range.end: {coverage_p_lin.y_range.end}" )
            except:
                pass

            # Background MASK coordinates
            if len(bg_mask):
                if taxid in MASK_CACHE:
                    mask_source.data = MASK_CACHE[taxid]
                else:
                    offset=1
                    m_ref_list=[]
                    m_str_list=[]
                    m_end_list=[]
                    for ref in pu['ref'].unique().tolist():
                        (acc, leng, taxid, tag) = ref.split('|')
                        coverage_p_log.title.text = f"Loading masks for {ref}..."
                        coverage_p_lin.title.text = f"Loading masks for {ref}..."
                        if ref in bg_mask:
                            p = re.compile('1+')
                            bitstr = bin(bg_mask[ref]).replace('0b','')
                            offset += int(leng)-len(bitstr)
                            iterator = p.finditer(bitstr)
                            mask_rgns = [match.span() for match in iterator]
                            for rgn in mask_rgns:
                                m_ref_list.append(ref)
                                m_str_list.append(rgn[0]+offset)
                                m_end_list.append(rgn[1]+offset-1)
                        offset += int(leng)
                    
                    mask_source.data = dict(
                        m_ref = m_ref_list,
                        m_str = m_str_list,
                        m_end = m_end_list
                    )

                    # save to cache
                    MASK_CACHE[taxid] = mask_source.data

                # calculate the length of mask
                for rgn in zip(mask_source.data['m_str'], mask_source.data['m_end']):
                    mask_bp += rgn[1]-rgn[0]
                
                logconfig.bokeh_logger.debug( "[DEBUG] [MASK] %s -- masked region: %sbp."%(name, mask_bp) )

            # BSAT coordinates
            if taxid in BSAT_CACHE:
                bsat_source.data = BSAT_CACHE[taxid]
            else:
                offset=1
                note_color={}
                note_color['common to other genomes'] = C_SET1[2]
                note_color['AMR'] = C_SET1[1]
                note_color['Biotoxin'] = C_SET1[0]
                note_color['Virulence'] = C_SET1[3]
                bsat = dict( b_ref = [], b_str = [], b_end = [], b_col = [], b_not = [])

                for ref in pu['ref'].unique().tolist():
                    logconfig.bokeh_logger.debug( f"[DEBUG] [BSAT] checking BSAT coordinates for {ref}..." )
                    (acc, leng, taxid, tag) = ref.split('|')
                
                    # find available bsat file
                    bsat_fn = f"{r_path}/{acc}.bsat.bed"
                    if not os.path.isfile(bsat_fn):
                        bsat_fn = f"{db_path}/BSAT_markers/{acc}.bsat.bed"

                    if os.path.isfile(bsat_fn) and os.path.getsize(bsat_fn) > 0:
                        coverage_p_log.title.text = f"Loading BSAT coordinates for {ref}..."
                        coverage_p_lin.title.text = f"Loading BSAT coordinates for {ref}..."
                        #loading coverage files
                        pu = pd.read_csv( bsat_fn, sep='\t', header=None, names=['ref', 'str', 'end', 'note'] )
                        #pu = pu.drop_duplicates(subset=['ref', 'str', 'end'])
                        bsat['b_ref'] += [ref for x in pu['note'].tolist()]
                        bsat['b_str'] += (pu['str']+offset).tolist()
                        bsat['b_end'] += (pu['end']+offset).tolist()
                        bsat['b_col'] += [note_color[x] for x in pu['note'].tolist()]
                        bsat['b_not'] += pu['note'].tolist()
                        logconfig.bokeh_logger.debug( f"[DEBUG] [BSAT] {bsat_fn} loaded." )
                    else:
                        logconfig.bokeh_logger.debug( f"[DEBUG] [BSAT] {bsat_fn} not available." )

                    offset += int(leng)

                bsat_source.data = bsat
                BSAT_CACHE[taxid] = bsat
            
            # change coverage plot title
            coverage_p_log.title.text = f"{name} (mask: {mask_bp}bp)" if mask_bp>0 else name
            coverage_p_lin.title.text = f"{name} (mask: {mask_bp}bp)" if mask_bp>0 else name
        else:
            reset_genome_plot()
    else:
        reset_genome_plot()

def update_taxainfo(attr, old, new):
    idx = 0
    fields = ['taxa','score','read_rsnb','lnr_cov','doc','pathogen','p_loc','p_src','p_host','p_dse']
    values = ["--" for field in fields]

    if len(dp_source.data['taxa']):
        try:
            idx = new[0]
        except:
            pass

        values = []
        for field in fields:
            val = dp_source.data[field][idx] if dp_source.data[field][idx] else "--"
            values.append(val)

    global taxa_rec_div
    taxa_rec_div.text="""
    <div class="table">
    <div class="row header green">
        <div class="cell">NAME</div>
        <div class="cell">SCORE</div>
        <div class="cell">READ COUNT (rsnb)</div>
        <div class="cell">GENOME COVERAGE</div>
        <div class="cell">DEPTH OF COVERAGE</div>
    </div>
    <div class="row">
        <div class="cell" id='rec_taxa'>%s</div>
        <div class="cell" id='rec_score'>%s</div>
        <div class="cell" id='rec_read_rsnb'>%s</div>
        <div class="cell" id='rec_lnr_cov'>%s</div>
        <div class="cell" id='rec_doc'>%s</div>
    </div>
    <div class="row header">
        <div class="cell">PATHOGEN</div>
        <div class="cell">LOCATION</div>
        <div class="cell">SOURCE</div>
        <div class="cell">HOST</div>
        <div class="cell">DISEASE</div>
    </div>    
    <div class="row">
        <div class="cell" id='rec_pathogen'>%s</div>
        <div class="cell" id='rec_p_loc'>%s</div>
        <div class="cell" id='rec_p_src'>%s</div>
        <div class="cell" id='rec_p_host'>%s</div>
        <div class="cell" id='rec_p_dse'>%s</div>
    </div>
    </div>"""%tuple(values)

###############################################################################
# control panel
###############################################################################

cols = [
    'None',
    'Read count (raw)',
    'Read count (raw, only primaly align)',
    'Read count (RNR)',
    'Read count (RSNB)',
    'Genome coverage',
    'RPKM',
    'Depth of coverage (RSNR)',
    'Depth of coverage',
    'Score',
    'Score (uniq)',
    'Score (background)',
    'Relative abundance',
]

cols_y = [
    'Read count (raw)',
    'Read count (raw, only primaly align)',
    'Read count (RNR)',
    'Read count (RSNB)',
    'RPKM',
    'Depth of coverage (RSNR)',
    'Depth of coverage',
    'Relative abundance',
]

cols_s = [
    'None',
    'Genome coverage',
    'Score',
    'Depth of coverage (RSNR)',
    'Depth of coverage',
    'Relative abundance'
]

cols_p = [
    'None',
    'Genome coverage',
    'Score',
    'Relative abundance'
]

patho      = RadioButtonGroup(labels=["Pathogen only", "All taxonomies"], active=def_val_patho)
rank       = Select(title='Rank', value='species', options=ranks)
coff_btn   = Button(label="No / Default extra cutoffs", button_type="success")
#min_len    = Slider(title="Minimum linear length", start=0, end=500, value=def_val_min_len, step=1)
min_cov    = Slider(title="Minimum genome coverage", start=0, end=1, value=def_val_min_cov, step=0.01)
max_r_raw  = Slider(title="Minimum reads", start=0, end=500, value=def_val_max_r_raw, step=1)
max_r_rsnb = Slider(title="Minimum reads (rsnb)", start=0, end=100, value=def_val_max_r_rsnb, step=1)
min_score  = Slider(title="Minimum score", start=0, end=1, value=def_val_min_score, step=0.1)
min_dc     = Slider(title="Minimum depth coverage (mX)", start=0, end=10000, value=def_val_min_dc, step=10)
min_rsdc   = Slider(title="Minimum depth coverage (RS) (mX)", start=0, end=1000, value=def_val_min_rsdc, step=1)
add_cut_c  = Select(title="Additional cutoff field", value='None', options=cols)
add_cut_v  = TextInput(title="Additional cutoff value")
org_n      = TextInput(title="Organism name contains")
dse_n      = TextInput(title="Disease name contains")
y          = Select(title='Y-Axis', value='Read count (RNR)', options=cols_y)
size       = Select(title='Size', value='Relative abundance', options=cols_s)
color      = Select(title='Color', value='Score', options=cols_p)
dl_btn     = Button(label="Export to CSV", button_type="success")

# control panel
controls = [patho, rank, coff_btn, min_cov, max_r_raw, max_r_rsnb, min_rsdc, min_dc, min_score, add_cut_c, add_cut_v,
            org_n, dse_n, y, size, color, dl_btn]

patho.on_change('active', lambda attr, old, new: update())

for control in [controls[1]]+controls[3:-1]:
    control.on_change('value', lambda attr, old, new: update())

dl_btn.callback = CustomJS(args=dict(source=dp_source), code=open(join(dirname(__file__), "download.js")).read())
coff_btn.on_click( cutoff_toggle )

ctrl_panel = widgetbox(*controls, sizing_mode='fixed', width=300)

# update taxa info when the user adjusts filters
dp_source.on_change('data', update_genome_cov )
dp_source.on_change('data', update_rankspec_rc )
dp_source.on_change('data', update_taxainfo)

# update taxa info when the user selects a taxa
dp_source.selected.on_change('indices', update_genome_cov )
dp_source.selected.on_change('indices', update_rankspec_rc )
dp_source.selected.on_change('indices', update_taxainfo)

###############################################################################
# LAYOUT
###############################################################################

header         = Div(text="""
    <H1>PanGIA Bioinformatics</H1>""", width=800)

# for top or selected record
taxa_rec_header = Div(text="""
    <H3>Selected organism / Top pathogen record</H3>
    <p>The detail of selected record will be displayed in this section.</p>""", width=800)
taxa_rec_div    = Div(text="""
    <div class="table">
    </div>""", width=770)
table_header   = Div(text="""
    <H3>Overview %s</H3>
    <p>PanGIA bioinformatics results are listed in the table below.</p>"""%proj_name, width=800)
plot_header    = Div(text="""
    <H3>Exploratory</H3>
    <p>This section will display results in a dot plot associated with the control panel dynamically. The chart displays organisms vs. read count / genome coverage / score by default. Cutoffs can be adjusted using the control panel on the right. Organism name is displayed at x-axis, read count (rnr) at y-axis, genome coverage as circle size and score as color.
        <div style="text-align: right;">
            Color: <table style="border-spacing: 0px; border-width: 1px; font-size: 0.8em; display: inline-table;">
            <tr>
            <td height="7px"> 0.00 </td>
            <td height="7px" width="7px" style="background-color: #313695;border: #444444 1px solid;border-right: None;"></td>
            <td height="7px" width="7px" style="background-color: #4575b4;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #74add1;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #abd9e9;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #e0f3f8;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #ffffbf;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #fee090;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #fdae61;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #f46d43;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #d73027;border-top: #444444 1px solid;border-bottom: #444444 1px solid;"></td>
            <td height="7px" width="7px" style="background-color: #a50026;border: #444444 1px solid;border-left: None;"></td>
            <td height="7px"> 1.00 </td>
            </tr>
            </table>
        </div>
    </p>""", width=800)

coverage_header = Div(text="""
    <span style="display: block; font-size: 1.17em; margin-bottom: 0px; font-weight: bold;">Genome Coverage</span>"""
, width=800)

footer = Div(text="""
    <div>PanGIA-VIS version %s</div>"""%__version__, width=1015)

l = layout([
    [header],
    [pieInReadsFigure, pieFlagFigure, piePathoFigure],
    [table_header],
    [result_table],
    [plot_header],
    [p, ctrl_panel],
    [coverage_header],
    [coverage_p_tabs],
    [taxa_rec_header],
    [taxa_rec_div, [rs_p,rsi_p]],
    [footer]
], sizing_mode='fixed')

setupPieFigures( piePathoFigure )
setupPieFigures( pieInReadsFigure )
setupPieFigures( pieFlagFigure )
update() # initial load of the data
loading_log() # initiate dashboard
loadBgMask(bgfiles) # load bg JSON file if possible

curdoc().add_root(l)
curdoc().title = "PanGIA Bioinformatics"