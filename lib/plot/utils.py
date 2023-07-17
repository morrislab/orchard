###############################################################################################################
# utils.py
# 
# Contains the source code for reading/writing 'Neutree' namedtuples 
###############################################################################################################
import sys, os, urllib, base64
from pyvis.network import Network
import numpy as np
import pandas as pd
from io import BytesIO

# matplotlib imports
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# omics data imports
from omicsdata.tree.parents import parents_to_adj
from omicsdata.tree.adj import adj_to_anc
from omicsdata.npz.archive import Archive
from omicsdata.ssm.parse import load_params

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "orchard"))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "metrics"))

from constants import KEYS_ORCH_NPZ, KEYS_PARAMS
import neutree 

# CONSTANTS
color_palette =  ["#E53939", "#E2E539", "#45E539", "#39E5CA", "#3986E5", 
                  "#9039E5", "#E539A1", "#B45555", "#9FB455", "#60B455", 
                  "#55B3B4", "#5568B4", "#B455A8"] 

img_html_style = 'style="width:100%"'
tree_html_header = '<h1 style="text-align:center">Tree</h1>'
tree_card_html_style = 'style="float:left; width:50%; height:100vh"'
right_panel_html_header = '<h1 style="text-align:center">Data</h1>'
right_panel_div_html_style = 'style="flex-direction:column;float:left;width:50%;height:90vh;display:flex;overflow:auto"'

def create_panel_html(tree_html, right_panel_html):
    """Method for modifying the tree html outputted by Pyvis to create two panels

    Parameters
    -----------
    tree_html : str
        the html for the Pyvis output
    right_panel_html : list
        a list of html snippets to place in the right panel

    Returns
    ---------
    str
        the html for the complete html file
    """
    tree_html = tree_html.replace('<div class="card" style="width: 100%">', '<div><div class="card" %s">%s' % (tree_card_html_style, tree_html_header))
    tree_html = tree_html.replace('</body>\n</html>', '%s</div></div></body>\n</html>' % right_panel_html_header)
    tree_html = tree_html.replace('</body>\n</html>', '<div %s>%s</div></div></body>\n</html>' % (right_panel_div_html_style, ''.join(right_panel_html)))
    return tree_html

def convert_fig_to_html(fig):
    """ Convert Matplotlib figure 'fig' into a <img> tag for HTML use using base64 encoding. """
    canvas = FigureCanvas(fig)
    fig_bytes = BytesIO()
    canvas.print_png(fig_bytes)
    data = fig_bytes.getvalue()

    encoded = base64.b64encode(data).decode('utf-8')
    return '<div %s><img src="data:image/png;base64,%s"></div>' % (img_html_style,encoded)

def load_npz_data(input_fn, params_fn):
    """Function for loading data either from an npz file or a Neutree file
    
    Parameters
    -----------
    input_fn : str
        A path to either an omicsdata Archive or a Neutree file

    Returns
    --------
    ndarray
        a list of lists where each sub-list is a parents vector
    ndarray
        a list of 2D numpy arrays where each 2D numpy array is the 
        cellular prevalence matrix for a particular parents vector/tree
    list
        a list of lists where each sub-list are the 'id' values for the mutations in the same cluster/clone
    list
        a list of sample names 
    """
    data = None
    structs, Fs, clusters, samples, llhs = None, None, None, None, None
    exceptions = []

    # try openining as an Archive
    try:
        data = Archive(input_fn)
        structs = data.get(KEYS_ORCH_NPZ.struct_key)
        Fs = data.get(KEYS_ORCH_NPZ.F_key)
        clusters = data.get(KEYS_ORCH_NPZ.clusters_key)
        samples = data.get(KEYS_ORCH_NPZ.sampnames)
        llhs = data.get(KEYS_ORCH_NPZ.llh_key)
    except Exception as e:
        exceptions.append(e)
    
    # try to open as a Neutree if we obtained an exceptionf
    if len(exceptions) > 0:
        try:
            assert isinstance(params_fn, str), "A parameter file is required for plotting when using a Neutree file"
            params = load_params(params_fn)
            ntree = neutree.load(input_fn)
            structs = ntree.structs
            Fs = ntree.phis
            clusters = ntree.clusterings[0] # assume all clusterings are the same
            samples = params[KEYS_PARAMS.samples_key]
            llhs = ntree.logscores
        except Exception as e:
            exceptions.append(e)
            print("Multiple exceptions occurred on loading %s" % input_fn)
            for e in exceptions:
                print(e)
            exit()

    return structs, Fs, clusters, samples, llhs

def calc_subpopfreq(parents, F):
    """Caculates the subpopulation frequency (\eta) of each clone in each sample
    
    Parameters
    -----------
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree
    F : ndarray
        a 2D numpy arrays where the rows are clones/mutations and the columns are samples, and each entry (i,s) is the cellular prevalence of 
        clone/mutation i in sample s that adhere to a perfect phylogeny (i.e., tree-constrained cellular prevalence)

    Returns
    --------
    ndarray
        a 2D numpy array of subpopulations frequencies for each clone/mutation. The rows are clones/mutations and the columns are samples, 
        and the entry (i,s) is the percentage of cells in that sample that have the mutational profile i, in sample s
    """
    adj = parents_to_adj(parents)
    anc = adj_to_anc(adj)
    anc_inv = np.linalg.inv(anc)
    return np.dot(anc_inv, F)

def CP_fig(F, samples, clusters, F_hat=None, color_palette=color_palette):
    """
    Generates a cellular prevalence plot for each clone/mutation in F

    Parameters
    -----------
    F : ndarray
        a 2D numpy arrays where the rows are clones/mutations and the columns are samples, and each entry (i,s) is the cellular prevalence of 
        clone/mutation i in sample s that adhere to a perfect phylogeny (i.e., tree-constrained cellular prevalence)
    samples : list
        a list of sample names
    clusters : list
        a list of lists where each sub-list are the 'id' values for the mutations in the same cluster/clone
    F_hat : ndarray, optional
        a 2D numpy arrays where the rows are clones/mutations and the columns are samples, and each entry (i,s) is the data-implied cellular prevalence of 
        clone/mutation i in sample s
    color_palette : list, optional
        a list of hex values for colors to give the different subplots

    Returns
    --------
    object
        a matplotlib.figure class instance
    """

    # determine the size, number of colors, and title for the cellular prevalence plot
    N, M = F.shape
    width = M+4
    plot_title = ""
    colors = [] # color for each subplot
    ylabels = [] # ylabel for each subplot

    # if we're given the data-implied cellular prevalence matrix, then we'll compute the interleaved cellular prevalence
    if isinstance(F_hat, np.ndarray): 
        height = N*3
        data = []
        for i in range(N):
            data.append(F_hat[i])
            data.append(F[i])
            data.append(np.abs(F[i]-F_hat[i]))
            colors.append(color_palette[i%len(color_palette)])
            colors.append(color_palette[i%len(color_palette)])
            colors.append("black")
        ylabels = ["data", "tree", "error"]*N
        plot_title="Interleaved cellular prevalence \n(Total Error = %.1f)" % np.abs(F-F_hat).sum()

    else: # otherwise, we'll just show the tree-constrained cellular prevalence
        height = N
        data = F
        ylabels = ["tree"]*N
        plot_title="Cellular prevalence"
        for i in range(height):
            colors.append(color_palette[i%len(color_palette)])
            
    fig = plt.figure(figsize=(width,height))
    gs = GridSpec(height,2, figure=fig, width_ratios=[1,min((M,10))])
    
    # add cellular prevalence subplots
    for i in range(height):
        ax = fig.add_subplot(gs[i,1])
        ax.bar(samples, data[i]*100, color=colors[i], width=1.0)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        ax.set_title("")
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel(ylabels[i])
        ax.set_ylim(0,105)
        ax.margins(x=0, tight=True)

        if i == 0:
            ax.set_title(plot_title, pad=20, wrap=True)
        if i < height - 1:
            ax.set_xticks("")
        else: 
            ax.tick_params(axis='x', labelrotation = 90)

        for p in ax.patches:
            ax.annotate("%.1f%%" % p.get_height(), (p.get_x() + 0.2, p.get_height()+3))


    # add cluster members subplots
    cluster_height = int(height / len(clusters))
    for i,C in enumerate([[]] + clusters):
        idx = i*cluster_height 
        ax = fig.add_subplot(gs[idx:idx+cluster_height,0])
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 0:
            ax.set_title("Mutations in node",pad=20)

        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_ylabel("Node %d" % i)
        ax.text(0.5,0.5,",".join(C),horizontalalignment='center', verticalalignment='center')

        # add horizontal lines between subplots
        trans = transforms.blended_transform_factory(fig.transFigure, ax.transAxes)
        line = Line2D([0.025,.09], [0,0], color='black', transform=trans, linewidth=0.9)
        fig.lines.append(line)

    gs.tight_layout(fig)

    return fig

def SP_fig(parents, F, samples, color_palette=color_palette, width=0.9):
    """Bar graph of subpopulation frequency per sample
    
    Parameters
    -----------
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree
    F : ndarray
        a 2D numpy arrays where the rows are clones/mutations and the columns are samples, and each entry (i,s) is the cellular prevalence of 
        clone/mutation i in sample s that adhere to a perfect phylogeny (i.e., tree-constrained cellular prevalence)
    samples : list
        a list of sample names
    color_palette : list, optional
        a list of hex values for colors to give the different subplots
    width : float, optional
        the width of each bar in the plot, must be a value in the range (0,1)

    Returns
    --------
    object
        a matplotlib.figure class instance
    """
    # calculates subpopulation frequency using cellular prevalence matrix + tree
    subpopfreq = calc_subpopfreq(parents, F)*100 
    subpopfreq_dict = {}
    N,M = F.shape
        
    # load in the subpopulation frequency data and color for each node into a dictionary
    for i in range(0,N):
        cl = color_palette[i%len(color_palette)]
        subpopfreq_dict["Node %d" % i] = (subpopfreq[i], cl)

    fig, ax = plt.subplots(figsize=(M,int(N/2)), layout="constrained")
    bottom = np.zeros(M)

    # plot each node's data
    for node,(f,color) in subpopfreq_dict.items():
        p = ax.bar(samples, f, width, label=node, bottom=bottom, color=color)
        bottom += f

    # format the plot
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)   
    ax.tick_params(axis='x', labelrotation = 90)
    ax.set_title("Subpopulation frequency in each sample", wrap=True)
    ax.legend(bbox_to_anchor=(1, 1))
    ax.margins(x=0, tight=True)

    # annotate the bar plot with its subpopulation frequency if it's large enough to fit
    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        if height > 8:
            ax.text(x+width/2, 
                    y+height/2, 
                    "%.1f%%" % height, 
                    horizontalalignment='center', 
                    verticalalignment='center')
    #fig.tight_layout()
    return fig

def draw_tree(adj, clusters):
    """
    Uses Pyvis to a draw a tree and return its html

    Parameters
    -----------
    adj : ndarray
        the adjacency matrix of the tree
    clusters : list
        a list of lists where each sub-list are the 'id' values for the mutations in the same cluster/clone

    Returns
    --------
    str
        an html string of the Pyvis tree
    """
    net = Network(directed=True, layout=True)

    # add nodes
    net.add_node(0, label="0", font='30px arial black', title="root")

    for i in range(1,adj.shape[0]):
        i_title = str(",".join(clusters[i-1]))
        net.add_node(i, label=str(i), font='30px arial black', title=i_title)

    # use the adjacency matrix to determine which nodes have edges between them
    for (i,j) in zip(*np.where(adj == 1)):
        i,j = int(i), int(j)

        # add edge
        net.add_edge(i,j)

    return net.generate_html()

def make_html_text(string, title="Description"):
    """Makes a description to add to the output html file
    
    Parameters
    ------------
    description : str
        either a text description, or an existing text file that can be read for a description

    Returns
    --------
    str
        an html snippet containing the description
    """
    text = string
    if os.path.exists(string):
        with open(string) as f:
            text = " ".join(f.readlines())

    return '<div style="padding-left:10px;word-wrap:break-word"><b>%s</b><br>%s</div>' % (title, text)
        
