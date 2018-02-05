import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''
Functions specific to coverage: reading in the data and creating the coverage track.
'''

def get_cov_df(cov_file, transcript, bp_start, bp_end, sample_l):

    '''
    Get the coverage data for a transcript. If bp_start and bp_end are None, then coverage data for the whole transcript will
    be extracted. The coordinates in cov_df and exon_coords_df are 1-based.
    
    Parameters
    ----------
    cov_file: str
        path to file containing the coverage data.
    transcript: str
        Ensembl transcript ID.
    bp_start: int
        base pair start coordinate.
    bp_end: int
        base pair end coordinate.
    sample_l: list of strs
        list of sample IDs.
    
    Returns
    -------
    cov_df: DataFrame
        contains the coverage data.
    '''    

    print("Reading in coverage data for " + str(len(sample_l)) + " samples...")
    if bp_start != None and bp_end != None:
        if bp_start >= bp_end:
            print("WARNING: bp_start " + str(bp_start) + " is not less than " + str(bp_end) + ".") 
    
    cov_df_chunker = pd.read_csv(cov_file, chunksize=1000)
    col_to_keep_l = ["chromStart","strand","position","cov","exon"]
    cov_df_chunk_l = []
    for cov_df_chunk in cov_df_chunker:
        cov_df_chunk = cov_df_chunk[cov_df_chunk["name"].str.contains(transcript)]
        cov_df_chunk["exon"] = cov_df_chunk["name"].str.split(pat=":").str.get(2)
        cov_df_chunk["cov"] = cov_df_chunk[sample_l].mean(axis=1)
        cov_df_chunk = cov_df_chunk[col_to_keep_l]
        #print cov_df_chunk
        cov_df_chunk_l.append(cov_df_chunk)
    
    cov_df = pd.concat(cov_df_chunk_l,ignore_index=True)
    del cov_df_chunk_l
    
    #Work out the bp and check whether bp_start <= bp <= bp_end
    #BED start positions are 0-based i.e. the first base in an exon is interpreted as +1 from the start position.
    #The end position is 1-based.
    cov_df["bp"] = cov_df["chromStart"] + cov_df["position"]
    cov_df.sort_values(by="bp", inplace=True)
    cov_df.index = range(len(cov_df.index))
    gene_strand = cov_df.iloc[0]["strand"]
    #cov_df.drop("strand", axis=1, inplace=True)
    cov_df["tp"] = range(1,len(cov_df.index)+1)
    if gene_strand == "-":
        cov_df["tp"] = range(len(cov_df.index),0,-1)
    cov_df.drop(["chromStart","position"], axis=1, inplace=True)
    
    return cov_df


def get_exon_coord_df(cov_df):

    '''
    Make an exon coordinate file.
    
    Parameters
    ----------
    cov_df: DataFrame
        contains the coverage information.
        
    Returns
    -------
    exon_coord_df: DataFrame
        contains the exon base pair and transcript position coordinates.
    '''

    #Exons: the exon start and end positions in this list are 1-based. 
    #cov_df["exon"] = cov_df["name"].str.split(pat=":").str.get(2)
    
    exon_coord_df = cov_df.groupby("exon")["bp"].agg([np.min,np.max])
    exon_coord_df.columns = ["start_bp","end_bp"] 
    
    #Add the transcript positions to exon_coords_df 
    get_tp_from_cov_df = lambda bp, cov_df: cov_df[cov_df["bp"] == bp].iloc[0]["tp"] 
    exon_coord_df["start_tp"] =  exon_coord_df["start_bp"].apply(func=get_tp_from_cov_df, cov_df=cov_df)
    exon_coord_df["end_tp"] =  exon_coord_df["end_bp"].apply(func=get_tp_from_cov_df, cov_df=cov_df)
    exon_coord_df.sort_values(by="start_tp", inplace=True)
    
    gene_strand = cov_df.iloc[0]["strand"]
    if gene_strand == "-":
        exon_coord_df.rename(columns={"start_bp":"end_bp", "end_bp":"start_bp", "start_tp":"end_tp", "end_tp":"start_tp"}, inplace=True)
    
    return exon_coord_df


def make_track(track, cov_df, bound_l, color_l, edge_color_l, setting_dict):
    
    '''
    Make the coverage track.
    
    Parameters
    ----------
    track: matplotlib.axes.Axes
        the axis for this coverage track.
    cov_df: DataFrame
        contains the coverage data.
    bound_l: list of ints
        bounds of the utrs and exons.
    color_l: list of strs
        colors for the utrs and exons.
    edge_color_l: list of strs
        edge colors for the utrs and exons.
    setting_dict: dictionary
        settings for making the png.'''
    
    for i in range(1,len(bound_l)):
        cov_in_bounds_df =  cov_df[(cov_df["tp"] >= bound_l[i-1]) & (cov_df["tp"] < bound_l[i])]
        #print cov_in_bounds_df
        track.fill_between(cov_in_bounds_df["tp"].tolist(), cov_in_bounds_df["cov"], facecolor=color_l[i-1], edgecolor=edge_color_l[i-1])

    track.set_xlabel('Position')
    plt.rc('text',usetex=True)
    track.set_ylabel(setting_dict["c_track_y_axis_label"], rotation="horizontal", size=setting_dict["c_track_fontsize"], ha='right', va='center')
    track.set_xlim((0,bound_l[-1]))
    track.set_ylim([0,cov_df["cov"].max()+10])
    track.grid(True)
    plt.setp(track.xaxis.get_ticklabels(), size=setting_dict["c_track_fontsize"])
    plt.setp(track.yaxis.get_ticklabels(), size=setting_dict["c_track_fontsize"])
    plt.setp(track.xaxis.get_label(), size=setting_dict["c_track_fontsize"])
    
    #NOTE, Above the xlim start is set to x[0] i.e. 1 so the track starts from 1, but the first xtick label in the plot is 0!
    #I have tried without success to replace this 0 with a 1.  
