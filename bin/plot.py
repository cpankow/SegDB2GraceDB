import matplotlib
matplotlib.use("agg")
from gwpy.plotter import SegmentPlot
from gwpy.segments import DataQualityDict, DataQualityFlag

def convert(known, active):
    new = DataQualityDict()
    for key in sorted(known):
        new[key] = DataQualityFlag(name=key, active=active[key], known=known[key])
    return new

def plot_seglist(known, active, gpstime):
    st, en = known.extent_all()
    segdict = convert(known, active)

    matplotlib.rcParams['ytick.labelsize'] = 8
    plot = SegmentPlot()
    ax = plot.gca()
    ax.plot(segdict)
    ax.set_epoch(int(round(gpstime, 6)))
    ax.set_xlim(st, en)
    ax.set_insetlabels(True)
    plot.save("test.png")
