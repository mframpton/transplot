import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
#import random
import numpy as np


'''
Functions specific to protein domains: reading in the protein domains and creating a protein domains track.
'''

def get_protein_domain_df(protein_domain_file, transcript_l, database, sortby_col_l):

    '''
    Read the protein domain information from a tsv file into a DataFrame.
    
    Parameters
    ----------
    protein_domain_file: str
        path to protein domain file.
    transcript_l: list of strs
        list of Ensembl transcript IDs.
    database: str
        protein domain database.
    sortby_col_l: list of strs
        columns to sort the DataFrame by
    
    Returns
    -------
    protein_domain_df: DataFrame
        contains the protein domain information.
    '''

    protein_domain_df = pd.DataFrame([])
    protein_domain_df_chunker = pd.read_csv(protein_domain_file, sep="\t", chunksize=1000)
    for protein_domain_df_chunk in protein_domain_df_chunker:
        protein_domain_df_chunk = protein_domain_df_chunk[["TranscriptID","Length","DomainID","Start","End","Domain_type","name","UniprotID"]]
        protein_domain_df_chunk = protein_domain_df_chunk[(protein_domain_df_chunk["TranscriptID"].isin(transcript_l)) & (protein_domain_df_chunk["Domain_type"]==database)]
        protein_domain_df = pd.concat([protein_domain_df, protein_domain_df_chunk], ignore_index=True)
    
    protein_domain_df.sort(sortby_col_l, inplace=True)

    return protein_domain_df


def get_protein_domain_color_s(protein_domain_df):
    
    '''
    Create a protein domain color series where the indexes are the protein domains, and the values are the colors.
    
    Parameters
    ----------
    protein_domain_df: DataFrame
        contains the protein domain information
    
    Returns
    -------
    protein_domain_color_s: Series
        contains colors indexed by protein domain.
    '''
    
    protein_domain_l = pd.unique(protein_domain_df["DomainID"])
    #random.shuffle(protein_domain_l)
    cmap = mpl.cm.get_cmap(name="Set1", lut=len(protein_domain_l))
    cmap_tmp_l = [cmap(i) for i in xrange(len(protein_domain_l))] 
    cmap_l = []
    for c in cmap_tmp_l:
        if isinstance(c, str):
            cmap_l.append(c)
        elif isinstance(c, tuple):
            cmap_l.append(",".join([str(f) for f in list(c)]))
    protein_domain_color_s = pd.Series(data=cmap_l, index=protein_domain_l)
    protein_domain_color_s.ix["none"] = "#D1D0CE"
    protein_domain_color_s.ix["utr"] = "white"
    
    return protein_domain_color_s


def make_track(track, protein_domain_df, utr_df, protein_domain_color_s, setting_dict, variant_track):

    '''
    Make the protein domain track.
    
    Parameters
    ----------
    track: matplotlib.axes.Axes
        axis for the protein domain track.
    protein_domain_df: DataFrame
        contains the protein domain information.
    utr_df: DataFrame
        contains the utr information.
    protein_domain_color_s: Series
        colors indexed by protein domain IDs.
    setting_dict: dictionary
        settings for making the png.
    variants_track: matplotlib.axes.Axes
        axis for the variant track.'''

    '''Work out the protein domain start and end transcript positions.'''
    protein_domain_df["start_tp_pc"] = (protein_domain_df["Start"]/protein_domain_df["Length"]).astype(float)
    protein_domain_df["end_tp_pc"] = (protein_domain_df["End"]/protein_domain_df["Length"]).astype(float)
    transcript_len = utr_df["end_tp"].max()
    tr_start = utr_df[utr_df["utr"]=="5'"]["end_tp"].max()
    tr_end = utr_df[utr_df["utr"]=="3'"]["start_tp"].min() - 1
    tr_diff = tr_end - tr_start
    protein_domain_df["start_tp"] = (tr_start+tr_diff*protein_domain_df["start_tp_pc"]).astype(int)
    protein_domain_df["end_tp"] = (tr_start+tr_diff*protein_domain_df["end_tp_pc"]).astype(int)
    
    '''Make the colorbar.'''
    [protein_domain_bound_l, bound_l, color_l] = get_bound_color_ls_for_cb(protein_domain_df, utr_df, protein_domain_color_s, 
                                                                           setting_dict["pd_track_stripe_min_bases"])
    cmap = mpl.colors.ListedColormap(color_l)
    norm = mpl.colors.BoundaryNorm(bound_l, cmap.N)
    cb = mpl.colorbar.ColorbarBase(track, cmap=cmap, norm=norm, boundaries=bound_l,
                                   spacing='proportional', orientation='horizontal', ticks=[], drawedges=False)#, alpha=0.5)
    plt.setp(track.get_yticklabels(),visible=False)
    plt.rc('text',usetex=True)
    track.set_ylabel(setting_dict["pd_track_y_label"], rotation='horizontal', ha='right', va=setting_dict["pd_track_y_label_va"],
                     position=(setting_dict["pd_track_y_label_x"], setting_dict["pd_track_y_label_y"]), size=setting_dict["pd_track_y_label_fontsize"])

    '''Make the protein domain track legend.'''
    make_legend(track, protein_domain_df, protein_domain_color_s, setting_dict["pd_track_legend_fontsize"], setting_dict["pd_track_legend_max_chars_per_row"],
                setting_dict["pd_track_legend_bbox"])
    
    '''Add tick marks for the protein domains to the variant track.'''
    if variant_track != None:
        variant_track = set_colorbar_ticks(variant_track, protein_domain_bound_l, transcript_len, "both", 5, 1)


def get_bound_color_ls_for_cb(protein_domain_df, utr_df, protein_domain_color_s, stripe_min_bases):
    
    '''
    Get the bound and color lists for the colorbar. These are (1) pd_bound_l, list of protein domain bounds
    which can be used to mark the protein domains on the variant track; (2) whole_trans_stripe_bound_l, list of
    regions covering the whole transcript in which contiguous regions have different contents (utr, protein
    domain(s), empty), and regions containing > 1 protein domain are striped. This list will be used to make the 
    protein domain colorbar; (3) whole_trans_stripe_color_l, list of colors corresponding to 
    whole_trans_stripe_bound_l.
    
    Parameters
    ----------
    protein_domain_df: DataFrame
        contains the protein domain information.
    utr_df: DataFrame
        contains the utr information.
    protein_domain_color_s: Series
        colors indexed by protein domain IDs.
    stripe_min_bases:
        minimum stripe width in bases.
    
    Returns
    -------
    pd_bound_l: list of ints
        list of protein domain bounds
    whole_trans_stripe_bound_l: list of ints
        list of utr and protein domain bounds
    whole_trans_stripe_color_l: list of strs and tuples of floats
        list of corresponding colors for the latter.
    '''
    
    '''Get the bounds without the utrs.'''
    protein_domain_start_tp_l = protein_domain_df["start_tp"].tolist()
    protein_domain_end_tp_l = protein_domain_df["end_tp"].tolist()
    pd_bound_l = [[protein_domain_start_tp_l[i], protein_domain_end_tp_l[i]] for i in xrange(len(protein_domain_start_tp_l))]
    pd_bound_l = [bound for bound_l in pd_bound_l for bound in bound_l]
    
    '''Add the utrs to the protein domains.'''
    bound_color_df = pd.concat([protein_domain_df, utr_df])
    bound_color_df = bound_color_df[["DomainID", "start_tp", "end_tp"]]
    bound_color_df["DomainID"].fillna(value="utr", inplace=True)
    bound_color_df.sort_values(["start_tp", "end_tp", "DomainID"], inplace=True)
    #print bound_color_df
    
    '''Divide the whole transcript into regions such that contiguous regions differ wrt their contents (utr, protein domain(s), empty)'''
    domain_l = bound_color_df["DomainID"].tolist()
    start_bound_l = bound_color_df["start_tp"].tolist()
    end_bound_l = bound_color_df["end_tp"].tolist()
    pd_utr_bound_ll = [[start_bound_l[i], end_bound_l[i]] for i in xrange(len(start_bound_l))]
    whole_trans_bound_ll = get_whole_trans_bound_ll(pd_utr_bound_ll)
    whole_trans_domain_ll = get_whole_trans_domain_ll(whole_trans_bound_ll, pd_utr_bound_ll, domain_l)
    #print whole_trans_domain_ll
    #Add the none regions.
    whole_trans_domain_ll = [["none"] if len(whole_trans_domain_ll[i]) == 0 else whole_trans_domain_ll[i] for i in xrange(len(whole_trans_domain_ll))]
    #print whole_trans_domain_ll
        
    '''Add bounds for stripes in the regions containing > 1 protein domain.'''
    whole_trans_stripe_bound_l, whole_trans_stripe_color_l = [],[]
    for i in xrange(len(whole_trans_bound_ll)):
        if len(whole_trans_domain_ll[i]) == 1:
            whole_trans_stripe_bound_l.append(whole_trans_bound_ll[i][0])
            whole_trans_stripe_color_l.append(protein_domain_color_s.ix[whole_trans_domain_ll[i][0]])
        else:
            [stripe_bound_l, stripe_domain_l] = generate_stripe_ls(whole_trans_bound_ll, whole_trans_domain_ll, i, stripe_min_bases)
            whole_trans_stripe_bound_l.extend(stripe_bound_l[:-1])#Don't need the region end bound (use next region start bound in its place).
            whole_trans_stripe_color_l.extend([protein_domain_color_s.ix[domain] for domain in stripe_domain_l])
    whole_trans_stripe_bound_l.append(whole_trans_bound_ll[-1][1])
    
    #print whole_trans_stripe_bound_l
    #print len(whole_trans_stripe_bound_l)
    #print whole_trans_stripe_color_l
    #print len(whole_trans_stripe_color_l)

    return [pd_bound_l, whole_trans_stripe_bound_l, whole_trans_stripe_color_l]


def get_whole_trans_bound_ll(pd_utr_bound_ll): 
    
    '''
    Divide the whole transcript into regions such that contiguous regions differ wrt their contents (utr, protein domain(s), empty), and return their bounds.
    
    Parameters
    ----------
    pd_utr_bound_ll: list of list of ints
        the start and end bounds of the utrs and protein domains.
    
    Returns
    -------
    whole_trans_bound_ll: list of list of ints
        lists containing start and end bound coordinates covering the whole transcript, where contiguous regions differ wrt their contents (utr, protein domain(s), neither). 
    '''
    
    #print "pd_utr_bound_ll: " + str(pd_utr_bound_ll)
    se_ll = [["s","e"] for bound_l in pd_utr_bound_ll]
    se_l = [bound_type for se_l in se_ll for bound_type in se_l]
    bound_l = [bound for bound_l in pd_utr_bound_ll for bound in bound_l]
    idx_l = sorted(range(len(bound_l)), key=lambda k: bound_l[k])
    bound_l = sorted(bound_l)
    se_l = [se_l[idx] for idx in idx_l]
    
    #print bound_l
    #print se_l
    
    whole_trans_bound_ll = []
    for i in xrange(len(bound_l)-1):
        bound = []
        if se_l[i] == "s":
            bound.append(bound_l[i])
        else:
            bound.append(bound_l[i] + 1)
        if se_l[i+1] == "s":
            bound.append(bound_l[i+1] - 1)
        else:
            bound.append(bound_l[i+1])
        if i > 0:
            if bound[0] <= whole_trans_bound_ll[-1][1]: 
                bound[0] = whole_trans_bound_ll[-1][1] + 1
        if bound[1] >= bound[0]:
            whole_trans_bound_ll.append(bound)

    return whole_trans_bound_ll
  

def get_whole_trans_domain_ll(whole_trans_bound_ll, pd_utr_bound_ll, domain_l):

    '''
    Get the utr or protein domain(s) in each region defined in whole_trans_bound_ll.
    
    Parameters
    ----------
    whole_trans_bound_ll: list of lists of ints
        lists containing start and end bound coordinates covering the whole transcript, where contiguous regions differ wrt their contents (utr, protein domain(s), empty). 
    pd_utr_bound_ll: list of lists of ints
        the start and end bounds of the utrs and protein domains.
    domain_l: list of strs
        list of protein domains corresponding to the bounds in pd_utr_bound_ll.
    
    Returns
    -------
    whole_trans_domain_ll: list of lists of strs
        the contents (utr, protein domain(s), empty) of each region in whole_trans_bound_ll.
    '''

    whole_trans_domain_ll = []
    for i in range(len(whole_trans_bound_ll)):
        overlap_domain_l = []
        new_set = set(xrange(whole_trans_bound_ll[i][0],whole_trans_bound_ll[i][1]+1))
        for j in range(len(pd_utr_bound_ll)):
            old_set = set(xrange(pd_utr_bound_ll[j][0],pd_utr_bound_ll[j][1]+1))
            if len(new_set.intersection(old_set)) > 0:
                overlap_domain_l.append(domain_l[j])
        whole_trans_domain_ll.append(overlap_domain_l)
    #print whole_trans_domain_ll
        
    return whole_trans_domain_ll


def generate_stripe_ls(whole_trans_bound_ll, whole_trans_domain_ll, i, stripe_min_bases):
    
    '''
    Generate the stripe lists for regions which contain greater than 1 protein domain. The 1st color will be for
    a protein domain not in the left-neighbouring region, and the last will be for a protein domain not in the
    right-neigbouring region.
    
    Parameters
    ----------
    whole_trans_bound_ll: list of lists of ints
        lists containing start and end bound coordinates covering the whole transcript, where contiguous regions differ wrt their contents (utr, protein domain(s), empty).
    whole_trans_domain_ll: list of lists of strs
        the contents (utr, protein domain(s), empty) of each region in whole_trans_bound_ll.
    i: int
        index in whole_trans_bound_ll.
    stripe_min_bases:
        minimum stripe width in bases.
    
    Returns
    -------
    stripe_bound_l: list of ints
        bounds for stripes, including end bound of region.
    stripe_domain_l: list of strs
        protein domains corresponding to the stripe bounds (alternates through the protein domains in this region).
    '''
    
    def get_domain_removed(domain_l_1, domain_l_2): 
        diff_l =  list(set(domain_l_1) - set(domain_l_2))
        new_domain = None if len(diff_l) == 0 else diff_l[0] 
        return new_domain
    
    '''Work out the 1st domain, last domain and other domains.'''
    first_domain = None if i == 0 else get_domain_removed(whole_trans_domain_ll[i], whole_trans_domain_ll[i-1])
    last_domain = None if i >= len(whole_trans_domain_ll) - 1 else get_domain_removed(whole_trans_domain_ll[i], whole_trans_domain_ll[i+1])
    other_domain_l = [domain for domain in whole_trans_domain_ll[i] if domain not in [first_domain, last_domain]]
    domain_l = filter(lambda x: x != None, pd.unique([first_domain] + other_domain_l + [last_domain])) 
    
    '''Work out the number of intervals.'''
    num_intervals = len(filter(lambda x: x != None, [first_domain] + other_domain_l + [last_domain]))
    interval_inc = num_intervals - 1 if first_domain == last_domain and first_domain != None else num_intervals
    region_len = whole_trans_bound_ll[i][1] - whole_trans_bound_ll[i][0]
    while float(region_len)/num_intervals >= stripe_min_bases:
        num_intervals += interval_inc
    
    '''Generate the stripe bounds and corresponding domains.'''
    stripe_bound_l = list(np.linspace(whole_trans_bound_ll[i][0], whole_trans_bound_ll[i][1], num_intervals+1))
    stripe_domain_l = [domain_l[k % len(domain_l)] for k in xrange(len(stripe_bound_l)-1)] 
    
    return [stripe_bound_l, stripe_domain_l]


def make_legend(protein_domain_track, protein_domain_df, protein_domain_color_s, fontsize, max_chars_per_row, bbox):
    
    '''
    Make the protein domains track legend.
    
    Parameters
    ----------
    protein_domain_track: matplotlib.axes.Axes
        axis for the protein domain track.
    protein_domain_df: DataFrame
        contains the protein domain information.
    protein_domain_color_s:
        colors indexed by protein domain IDs.
    fontsize: int
        fontsize for text in the legend.
    max_chars_per_row:
        maximum number of characters for 1 row of the legend.
    bbox:
        where to anchor the legend (uses the axes coordinate system so (0,0) is bottom left and (1.0,1.0) is top right.
    '''
    
    protein_domain_df = protein_domain_df.drop_duplicates(["DomainID"])
    pds_df_minus_utrs = protein_domain_df[~protein_domain_df["DomainID"].isin(["utr_3","utr_5"])]
    handles = [mpl.patches.Patch(facecolor=protein_domain_color_s.ix[protein_domain]) for protein_domain in pds_df_minus_utrs["DomainID"]]
    labels = [description for description in pds_df_minus_utrs["name"].dropna().tolist()]
    n_cols = get_legend_num_cols(pds_df_minus_utrs["name"].dropna().tolist(),max_chars_per_row)
    legend = protein_domain_track.legend(handles, labels, bbox_to_anchor=bbox, ncol=n_cols, mode="expand", bbox_transform=protein_domain_track.transAxes)
    
    for label in legend.get_texts():
        label.set_fontsize(fontsize)
        

def get_legend_num_cols(protein_domain_description_l, max_num_chars_per_row):

    '''
    Work out the number of columns to use in the protein domain track legend.
    
    Parameters
    ----------
    protein_domain_description_l: list of strs
        protein domain descriptions
    max_num_chars_per_row: int
        maximum number of characters for 1 row of the legend.
    
    Returns
    -------
    num_cols: int
        number of columns.
    '''

    patch_len = 6
    min_patch_text_gap = 2
    min_text_patch_gap = 5
    num_cols = len(protein_domain_description_l)
    while num_cols >= 1:
        zero_rows_too_long = True
        num_chars_in_row = 0
        for i in xrange(len(protein_domain_description_l)):
            num_chars_in_row += patch_len + min_patch_text_gap + len(protein_domain_description_l[i]) + min_text_patch_gap
            if num_chars_in_row > max_num_chars_per_row:
                zero_rows_too_long = False 
            if (i+1) % num_cols == 0:
                num_chars_in_row = 0 
        if zero_rows_too_long == True:
            break
        num_cols -= 1
        
    return num_cols


def set_colorbar_ticks(track, bound_l, transcript_len, top_or_bottom, markersize, markeredgewidth):
    
    '''
    Set the ticks on a colorbar.
    
    Parameters
    ----------
    track: matplotlib.axes.Axes
        axis for a track.
    bound_l: list of ints
    transcript_len: int
        length of the transcript.
    top_or_bottom: str
        indicates whether the ticks should be below or above or both.
    markersize: int
        tick size.
    markeredgewidth: int
        tick width.
    
    Returns
    -------
    track: matplotlib.axes.Axes
        the axis with ticks added.
    '''
    
    for i in xrange(len(bound_l)):
        bound_l[i] = float(bound_l[i]-1)/float(transcript_len-1) #This assumes the xlim is (1, transcript_length).
    
    plt.setp(track.xaxis.get_ticklines(),'markersize', markersize)
    plt.setp(track.xaxis.get_ticklines(),'markeredgewidth', markeredgewidth)
        
    track.xaxis.set_ticks(bound_l)
    track.xaxis.set_ticks_position(top_or_bottom)
    track.set_xticklabels([""]*len(bound_l))
    
    return track
