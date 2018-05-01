import utrs as u
import variants as v
from transplot import protdomains as pds
import coverage as c
import matplotlib.pyplot as plt
import sys
import pandas as pd


def make_protein_domain_color_file(protein_domain_file, transcript_l, database, sortby_col_l, out_path):

    '''
    Make a protein domain color file.
    
    protein_domain_file: str
        path to file containing the protein domain information.
    transcript_l: list of strs
        list of Ensembl transcript IDs.
    database: str
        protein domain database.
    sortby_col_l: list of strs
        list of names of columns to sort the DataFrame by.
    out_path: str
        path to write the protein domain color file to.
    '''
    
    print("make_protein_domain_color_file")
    protein_domain_df = pds.get_protein_domain_df(protein_domain_file, transcript_l, database, sortby_col_l)
    protein_domain_color_s = pds.get_protein_domain_color_s(protein_domain_df)
    protein_domain_color_s.to_csv(out_path)
    print("Written protein domain color file to {0}\n".format(out_path))
    

def make_exon_coord_file(cov_file, transcript, out_path):
    
    '''Make an exon coordinate file from a coverage file.
    
    cov_file: str
        path to coverage file.
    transcript: str
        Ensembl transcript ID.
    out_path: str
        path to write the exon coordinate file to.
    '''

    print("make_exon_coord_file")
    cov_df = c.get_cov_df(cov_file, transcript, None, None, [])
    exon_coord_df = c.get_exon_coord_df(cov_df)
    exon_coord_df.to_csv(out_path, index=True)
    print("Written exon coordinate file to {0}\n".format(out_path))
    

def make_png(transcript_l, title_l, track_l, sample_ll, utr_file_l, exon_coord_file_l, cov_file_l,
             variant_file_l, protein_domain_file_l, protein_domain_color_file, setting_dict, png_file):
    
    '''
    Make a png which contains coverage/variants/protein domain tracks for 1 or more transcripts, subject to space limitations.
    
    transcript_l: list of strs
        Ensembl transcript ID(s)
    title_l: list of strs
        title to use for each transcript.
    track_l: list of strs
        list of strings of length 3 which encode whether to generate each of the 3 tracks (coverage, variants, protein domains). 
    sample_ll: list of list of strs
        list of lists of sample IDs.
    utr_file_l: list of strs
        list of utr file paths.
    exon_coord_file_l: list of strs
        list of exon coordinate file paths.
    cov_file_l: list of strs
        list of coverage file paths.
    variant_file_l: list of strs
        list of variant file paths.
    protein_domain_file_l: list of strs
        list of protein domain file paths.
    protein_domain_color_file: str
        protein domain color file.
    setting_dict: dictionary
        settings for making the png.
    png_file: str
        path to write the png file to.
    '''

    print("make_png")

    '''Check the parameters are well-formed'''
    if len(transcript_l) == len(title_l) == len(track_l) == len(sample_ll) == len(utr_file_l) == len(exon_coord_file_l) == len(cov_file_l) == len(variant_file_l) == len(protein_domain_file_l):
        True
    else:
        print("ERROR: Parameters of make_png function which are lists must all be the same length.\n")
        return False
        
    '''Check that the maximum number of tracks is not exceeded.'''
    num_rows, previous_track = 0, None
    track_gap_dict = dict(zip(["c","v","pd", None],
                              [setting_dict["c_track_gap_rows"], setting_dict["v_track_gap_rows"], setting_dict["pd_track_gap_rows"], 0]))
    for track_s in track_l:
        if track_s[0] == "1":
            num_rows += track_gap_dict[previous_track] + setting_dict["c_track_rows"]
            previous_track = "c"
        if track_s[1] == "1":
            num_rows += track_gap_dict[previous_track] + setting_dict["v_anns_top_rows"] + setting_dict["v_track_rows"] + setting_dict["v_anns_bot_rows"] + setting_dict["v_key_rows"]
            previous_track = "v"
        if track_s[2] == "1":
            num_rows += setting_dict["pd_track_rows"] + track_gap_dict[previous_track]
            previous_track = "pd"
    if num_rows > setting_dict["fig_num_rows"]:
        print("ERROR: PNG figure requires {0} rows but there are only {1}\n".format(num_rows,setting_dict["fig_num_rows"]))
        return False

    '''Initialise the figure.'''
    fig = plt.figure()
    plt.rc('text', usetex=True)

    '''Make the axes.'''
    get_tp_from_exon_bp_l = lambda bp, exon_bp_l: exon_bp_l.index(int(bp)) + 1
    num_rows = setting_dict["fig_num_rows"]
    start_row = 0
    title_1_coords = None

    for i in range(len(transcript_l)):
        
        print("Transcript: {0}".format(transcript_l[i]))
        
        '''Read in the exon positions and the utrs.'''
        exon_coord_df = pd.read_csv(exon_coord_file_l[i], index_col="exon")
        strand = "+" if (exon_coord_df["end_bp"] > exon_coord_df["start_bp"]).all() else "-"
        if strand == "+":
            print("Transcription direction: forward")
        elif strand == "-":
            print("Transcript direction: reverse") 
        exon_bp_l = []
        for j in range(len(exon_coord_df.index)):
            if strand == "+":
                exon_bp_l.extend(range(exon_coord_df["start_bp"].iloc[j],exon_coord_df["end_bp"].iloc[j]+1))
            elif strand == "-":
                exon_bp_l.extend(range(exon_coord_df["start_bp"].iloc[j],exon_coord_df["end_bp"].iloc[j]-1,-1))
        utr_df = u.get_utr_df(utr_file_l[i], strand, transcript_l[i])
        #Add the transcript positions to utr_df. 
        utr_df["start_tp"] = utr_df["start_bp"].apply(func=get_tp_from_exon_bp_l, exon_bp_l=exon_bp_l)
        utr_df["end_tp"] = utr_df["end_bp"].apply(func=get_tp_from_exon_bp_l, exon_bp_l=exon_bp_l)
        
        '''Make the title track'''
        if i == 0:
            plt.figtext(setting_dict["title_1_fig_x"], setting_dict["title_1_fig_y"], title_l[i], fontsize=setting_dict["title_fontsize"])
            title_1_coords = fig.transFigure.transform((setting_dict["title_1_fig_x"], setting_dict["title_1_fig_y"]))
        else:
            title_track = plt.subplot2grid((num_rows,1), (start_row,0))
            title_track.set_axis_off()
            inv = title_track.transData.inverted()
            title_track.text(inv.transform(title_1_coords)[0], setting_dict["title_2_ax_y"], title_l[i], fontsize=setting_dict["title_fontsize"])
            start_row += setting_dict["t_track_rows"]
        
        '''Make the coverage track.'''
        if track_l[i][0] == "1":
            print("Making coverage track.")
            cov_df = c.get_cov_df(cov_file_l[i], transcript_l[i], None, None, sample_ll[i])
            [bound_l, color_l, edge_color_l] = get_exon_bound_color_l(exon_coord_df, utr_df, strand)
            coverage_track = plt.subplot2grid((num_rows,1), (start_row, 0), rowspan=setting_dict["c_track_rows"])
            start_row += setting_dict["c_track_rows"]
            c.make_track(coverage_track, cov_df, bound_l, color_l, edge_color_l, setting_dict)
            start_row += setting_dict["c_track_gap_rows"]
        
        '''Make the variants track.'''
        variant_track = None
        if track_l[i][1] == "1": #Make the variants track.
            print("Making variant track.")
            variant_df = v.get_variant_df(transcript_l[i], variant_file_l[i])
            variant_df.rename(columns={"pos":"bp"}, inplace=True)
            variant_df["tp"] = variant_df["bp"].apply(func=get_tp_from_exon_bp_l, exon_bp_l=exon_bp_l)
            #variant_df.drop_duplicates(subset=["GENE_prot_change","GENE_DNA_change"], inplace=True)
            variant_df.sort_values(by="tp", inplace=True)
            start_row += setting_dict["v_anns_top_rows"]
            variant_track = plt.subplot2grid((num_rows,1), (start_row,0))
            start_row += setting_dict["v_track_rows"]
            start_row += setting_dict["v_anns_bot_rows"]
            variant_key = plt.subplot2grid((num_rows,1), (start_row,0), rowspan=setting_dict["v_key_rows"])
            start_row += setting_dict["v_key_rows"]
            [bound_l, color_l, edge_color_l] = get_exon_bound_color_l(exon_coord_df, utr_df, strand)
            v.make_track(variant_track, bound_l[-1], bound_l, color_l, edge_color_l, variant_df, setting_dict, variant_key)
            start_row += setting_dict["v_track_gap_rows"]
        
        '''Make the protein domain track.'''
        if track_l[i][2] == "1": #Make the protein domains track.
            print("Making protein domain track.")
            protein_domain_df = pds.get_protein_domain_df(protein_domain_file_l[i], [transcript_l[i]], "Pfam", ["Start","End"])
            protein_domain_track = plt.subplot2grid((num_rows,1), (start_row,0))
            start_row += setting_dict["pd_track_rows"]
            protein_domain_color_df = pd.read_csv(protein_domain_color_file, header=None, names=["Domain", "Color"], index_col="Domain")
            protein_domain_color_s = pd.Series(data=protein_domain_color_df["Color"], index=protein_domain_color_df.index)
            del protein_domain_color_df
            protein_domain_color_s = protein_domain_color_s.apply(lambda x: x if "," not in x else tuple([float(f) for f in x.split(",")]))
            pds.make_track(protein_domain_track, protein_domain_df, utr_df, protein_domain_color_s, setting_dict, variant_track)
            start_row += setting_dict["pd_track_gap_rows"]

    fig.set_size_inches(setting_dict["fig_width_inches"], setting_dict["fig_height_inches"])
    plt.savefig(png_file, dpi=setting_dict["fig_dpi"])
    print("Written {0}.\n".format(png_file)) 


def get_exon_bound_color_l(exon_coord_df, utr_df, strand):

    '''
    Get the bounds, colors and edge colors required to generate a color bar that displays the utrs and exons for a transcript. 
    
    Parameters
    ----------
    exon_coord_df: DataFrame
        contains the exon base pair and transcript position coordinates.
    utr_df: DataFrame
        contains the utr base pair and transcript position coordinates.
    strand: str
        specifies whether the transcript is on the positive or negative strand.
    
    Returns
    -------
    exon_bound_color_ll: list of lists of ints and strs
        contains bound_l, the list of bounds, color_l, the list of colors and edge_color_l, the list of edge colors.
    '''

    bound_color_df = pd.concat([exon_coord_df,utr_df])
    bound_color_df.drop(["start_bp","end_bp"], axis=1, inplace=True)
    bound_color_df["utr"].fillna(value="exon", inplace=True)
    bound_color_df.sort_values(by=["start_tp","end_tp","utr"], inplace=True)
    bound_color_df.index = range(len(bound_color_df))
    bound_color_df.index.rename("",inplace=True)
    #print bound_color_df
    
    for i in range(1,len(bound_color_df)):
        if bound_color_df.loc[i,"end_tp"] == bound_color_df.loc[i-1,"end_tp"]:
            continue
        if bound_color_df.loc[i,"start_tp"] < bound_color_df.loc[i-1,"end_tp"]:
            bound_color_df.loc[i,"start_tp"] = bound_color_df.loc[i-1,"end_tp"] + 1    
    bound_color_df.drop_duplicates(["start_tp","end_tp"],inplace=True)
    #print bound_color_df
    
    bound_color_df["color"] = ["white"] * len(bound_color_df.index)
    exon_color_l = ["red","#6E6E6E"]
    exon_color_i = 0
    for i in bound_color_df.index.tolist():
        #print i
        if bound_color_df.loc[i, "utr"] == "5'" or bound_color_df.loc[i,"utr"] == "3'":
            bound_color_df.loc[i,"color"] = "white"
        else:
            bound_color_df.loc[i,"color"] = exon_color_l[exon_color_i]
            exon_color_i = 1 - exon_color_i
    
    bound_l = bound_color_df["start_tp"].tolist()
    transcript_len = bound_color_df.iloc[len(bound_color_df.index)-1]["end_tp"]
    if bound_l[-1] != transcript_len:
        bound_l.append(transcript_len)
    color_l = bound_color_df["color"].tolist()

    edge_color_l = ["black" if color == "white" else color for color in color_l]

    exon_bound_color_ll = [bound_l, color_l, edge_color_l]
    
    return exon_bound_color_ll

