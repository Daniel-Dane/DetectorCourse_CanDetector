from rootpy.plotting import Hist
from matplotlib.font_manager import FontProperties

font0 = FontProperties()
font = font0.copy()
font.set_style('italic')
font.set_weight('bold')
font.set_size('x-large')
font_wip = font0.copy()
font_wip.set_style('italic')
font_wip.set_weight('medium')

# load from mca file
def mca_to_hist(filename, do_print = True):
    roi_on = False
    data_on = False
    dpp_on = False
    data_line = 1
    nbins = 1024
    r_min = 0
    r_max = 0
    time = 0

    h = Hist(nbins, 0, nbins)

    import sys
    enc = {} if sys.version_info[0]==2 else {"encoding":"latin-1"}
    for line in open(filename, 'r', **enc):
        if roi_on:
            val = line.split()
            r_min = val[0]
            r_max = val[1]
            roi_on = False

        if data_on:
            if line.find("END") > -1:
                data_on = False
                continue
            h.SetBinContent(data_line, float(line.split()[0]))
            data_line += 1

        if dpp_on:
            if line.find("END") > -1:
                dpp_on = False
                continue
            if line.find("Real Time: ") > -1:
                time = float(line[11:])

        if line.find("ROI") > -1:
            roi_on = True

        if line.find("DATA") > -1:
            data_on = True

        if line.find("DPP STATUS") > -1:
            dpp_on = True

    if do_print:
        print("region of interest: {:} to {:}".format(r_min, r_max))

    return (h, r_min, r_max, time)

# show text on figure below title
def show_text( text, ax, x=0.05, y=0.9, verticalalignment='bottom', horizontalalignment='left', fontproperties=font_wip, ha="left" ) :
    ax.text(x, y, text, verticalalignment=verticalalignment, horizontalalignment=horizontalalignment, fontproperties=fontproperties, transform=ax.transAxes, ha=ha)

# show 1-Group title on figure
show_title_subtitle=""
def show_title( ax, x=0.05, y=0.9, verticalalignment='bottom', horizontalalignment='left', fontproperties=font ) :
    ax.text(x, y, '1-Group', verticalalignment=verticalalignment, horizontalalignment=horizontalalignment, fontproperties=fontproperties, transform=ax.transAxes)
    if show_title_subtitle:
        if x < 0.5:
            show_text( show_title_subtitle, ax, x+0.2, y, verticalalignment, horizontalalignment )
        else:
            show_text( show_title_subtitle, ax, x-0.3, y, verticalalignment, horizontalalignment )