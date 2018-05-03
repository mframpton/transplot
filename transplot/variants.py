import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import regex as re
import sys

'''
Functions specific to variants: reading in the variant information and creating a variants track.
'''

def get_variant_df(transcript, variant_file):
    
    '''Read the variant information from a tsv file into a DataFrame.
    
    Args:
        | transcript (str): Ensembl transcript ID.
        | variant_file (str): path to file containing the variants.
    
    Returns:
        variant_df (DataFrame): contains the variant information.
    '''
    
    variant_df = pd.DataFrame([])
    variant_df_chunker = pd.read_csv(variant_file, sep="\t", chunksize=1000)
    for variant_df_chunk in variant_df_chunker:
        variant_df_chunk = variant_df_chunk[variant_df_chunk["featureID"] == transcript]
        variant_df = pd.concat([variant_df, variant_df_chunk], ignore_index=True)
    
    variant_df.fillna("NULL", inplace=True)
    
    return variant_df


def make_track(variant_track, transcript_len, bound_l, color_l, edge_color_l, variant_df, setting_dict, variant_key):

    '''Make the variants track.
    
    Args:
        | variant_track (matplotlib.axes.Axes): axis for the variants track
        | exon_coord_df (DataFrame): contains the exon base pair and transcript position coordinates.
        | bound_l (list of ints): contains the exons and utr bounds to be used in making the color bar.
        | color_l (list of strs): contains the list of colors to be used in making the color bar.
        | edge_color_l (list): contains the list of edge colors to be used in making the color bar.
        | variant_df (DataFrame): contains the variants information.
        | setting_dict (dictionary): settings for making the png.
        | variant_key (matplotlib.axes.Axes): axis for the variants key.
    '''

    #(1) Make the color bar.
    cmap = mpl.colors.ListedColormap(color_l)
    norm = mpl.colors.BoundaryNorm(bound_l, cmap.N)
    cb = mpl.colorbar.ColorbarBase(variant_track, cmap=cmap, norm=norm, boundaries=bound_l, spacing='proportional',
                                   orientation='horizontal', drawedges=False)
    plt.setp(variant_track.get_xticklabels(), visible=False)
    plt.rc('text', usetex=True)
    variant_track.set_ylabel(setting_dict["v_track_y_axis_label"], rotation='horizontal', ha='right', va='center', size=setting_dict["v_track_fontsize"])
    
    #(2) Annotate variant track with variants.
    annotate_track_with_variants(variant_track, variant_df, transcript_len, setting_dict)
    
    #(3) Make the variant annotations key.
    make_variant_annotations_key(variant_key, variant_df, setting_dict)


def annotate_track_with_variants(variant_track, variant_df, transcript_len, setting_dict):

    '''Annotate the variants track with arrows for the variants.
    
    Args:
        | variant_track (matplotlib.axes.Axes): axis for the variant track.
        | variant_df (DataFrame): contains the variant information.
        | transcript_len (int): transcript length.
        | setting_dict (dictionary): settings for making the png.
        
    Returns:
        variant_track (matplotlib.axes.Axes): axis for the variant track.
    '''

    if variant_df.shape[0] == 0:
        return variant_track
    
    #Add columns to variant_df for annotating variants with arrows: ID, axes x coordinates, arrow bin.
    variant_df["id"] = pd.Series([str(i) for i in range(1,variant_df.shape[0]+1)], index=variant_df.index)
    variant_df["top"] = variant_df.apply(lambda x: 1 if setting_dict["v_track_vars_t_or_b"][x["effect"]] == "T" else 0,
                                         axis=1) #Column for whether variant should be annotated with a top or bottom arrow.
    #Add column to determine whether each arrow should be merged with the previous.
    variant_df["trans_pos_pc"] = variant_df["tp"]/transcript_len #axes x coordinates.
    get_display_from_axes_coords = lambda axes_x_coord: variant_track.transAxes.transform((axes_x_coord,0))[0] - variant_track.transAxes.transform((0,0))[0]
    variant_df["num_pixels_diff"] = variant_df["trans_pos_pc"].diff().apply(get_display_from_axes_coords)
    variant_df["merge_wt_prev"] = variant_df["num_pixels_diff"] <= setting_dict["v_track_merge_pixel_thresh"]
    variant_df["arrow_bin"] = [0]*len(variant_df.index)
    arrow_bin = 1
    variant_df.at[0,"arrow_bin"] = arrow_bin
    for i in range(1,len(variant_df.index)):
        if variant_df.iloc[i]["merge_wt_prev"] == False or variant_df.iloc[i]["top"] != variant_df.iloc[i-1]["top"]:
            arrow_bin += 1
        variant_df.at[i,"arrow_bin"] = arrow_bin
    
    #Create new dataframe where each row corresponds to 1 arrow.
    x_pos_s = variant_df.groupby("arrow_bin")["trans_pos_pc"].mean()
    text_s = variant_df.groupby("arrow_bin").apply(lambda x: ",".join(x["id"].tolist())) #Make the annotation text string for a variant arrow.
    top_s = variant_df.groupby("arrow_bin").apply(lambda x: x["top"].tolist()[0]) #Get whether the arrow is a top or bottom arrow.   
    heights_s = get_arrow_height_s(top_s, setting_dict)
    variant_annotation_df = pd.concat([x_pos_s, text_s, top_s, heights_s],axis=1)
    variant_annotation_df.rename(columns={"trans_pos_pc":"x", 0:"text", 1:"top", 2:"height"}, inplace=True)

    #Annotate the variants.
    variant_annotation_df.apply(axis=1, func=annotate_track_with_arrow, variant_track=variant_track, setting_dict=setting_dict)    

    #Add text to indicate which variants types are annotated above and below the colorbar.
    variant_track.text(setting_dict["v_track_vars_text_top_x"], setting_dict["v_track_vars_text_top_y"], setting_dict["v_track_vars_text_top"],
                        ha='center', va='bottom', size=setting_dict["v_track_fontsize"])
    variant_track.text(setting_dict["v_track_vars_text_bot_x"], setting_dict["v_track_vars_text_bot_y"], setting_dict["v_track_vars_text_bot"],
                        ha='center', va='top', size=setting_dict["v_track_fontsize"])
    
    return variant_track


def get_arrow_height_s(top_s, setting_dict):
    
    '''Get the height of the arrow.
    
    Args:
        top_s (Series): indicates whether a variant is annotated with an arrow above or beneath the colorbar.
    
    Returns:
        height_s (Series): indicates the height the arrow used to annotate each variant.
    '''
    
    top_l = top_s.tolist()
    height_l = []
    t_bin,b_bin = 0,0
    for i in range(len(top_l)):
        if top_l[i] == 0:
            height_l.append(t_bin % setting_dict["v_track_num_arrow_heights"])
            t_bin += 1
        else:
            height_l.append(b_bin % setting_dict["v_track_num_arrow_heights"])
            b_bin += 1
    
    height_s = pd.Series(height_l, index=top_s.index)
    
    return height_s

    
def annotate_track_with_arrow(arrow_bin, variant_track, setting_dict):

    '''Annotate the variants track with an arrow.
    
    Args:
        | arrow_bin (int):
        | variant_track (matplotlib.axes.Axes): axis for the variant track.
        | setting_dict (dictionary): settings for making the png.
    '''

    if arrow_bin["top"] == 1:
        variant_track.annotate(arrow_bin["text"], xy=(arrow_bin["x"], setting_dict["v_track_t_arrow_head"]), xycoords='data',
                                xytext=(arrow_bin["x"], setting_dict["v_track_t_start"]+arrow_bin["height"]*setting_dict["v_track_t_inc"]), textcoords='data',
                                ha='center', arrowprops=dict(arrowstyle=setting_dict["v_track_arrow_style"], alpha=0.75), fontsize=8)
    else:
        variant_track.annotate(arrow_bin["text"], xy=(arrow_bin["x"], setting_dict["v_track_b_arrow_head"]), xycoords='data',
                                xytext=(arrow_bin["x"], setting_dict["v_track_b_start"]-arrow_bin["height"]*setting_dict["v_track_b_inc"]), textcoords='data',
                                ha='center', arrowprops=dict(arrowstyle=setting_dict["v_track_arrow_style"], alpha=0.75), fontsize=8)
    

def make_variant_annotations_key(variant_key, variant_df, setting_dict):
    
    '''Make a variant annotations key.
    
    Args:
        | variant_key (matplotlib.axes.Axes): axis for the variant key.
        | variant_df (DataFrame): contains the variant information.
        | setting_dict (dictionary): settings for making the png.
    
    Returns:
        variant_key (matplotlib.axes.Axes): axis for the variant key.
    '''
    
    plt.rc('text',usetex=True)
    variant_key.set_axis_off()
    variant_key_txt = r'''\begin{tabular}{''' + 'l'*setting_dict["v_key_num_cols"] + '''} \\\\ ''' 
    #variant_df["var_type_abbrev"] = variant_df.apply(lambda x: setting_dict["v_track_var_abbrevs"][x["effect"]] 
    #                                                 if setting_dict["v_track_var_abbrevs"].has_key(x["effect"]) != -1 else "", axis=1) #Abbreviate a variant type.
    variant_df["var_type_abbrev"] = variant_df.apply(lambda x: setting_dict["v_track_var_abbrevs"][x["effect"]] 
                                                     if "v_track_var_abbrevs" in setting_dict.keys() else "", axis=1) #Abbreviate a variant type.
    variant_df["variant_txt_str"] = variant_df["id"] + ": " + variant_df["dnachange"] + ", " + variant_df["prot_change"] + ", (" + variant_df["var_type_abbrev"] + ")"
    variant_df["variant_txt_str"] = variant_df["variant_txt_str"].map(mark_up_special_chars)    
    add_multicol = lambda cell_txt: "\multicolumn{2}{l}{" + cell_txt + "}" if len(cell_txt) > setting_dict["v_key_max_chars_per_col"] else cell_txt
    variant_df["variant_txt_str"] = variant_df["variant_txt_str"].map(add_multicol)
    var_txt_str_l = variant_df["variant_txt_str"].tolist()
    col_idx = 0
    for var_txt_str in var_txt_str_l:
        cols_for_str = 1
        m = re.search("\\multicolumn\{(\d*)\}",var_txt_str)
        if m != None: 
            cols_for_str = int(m.group(1))
        if col_idx + cols_for_str > setting_dict["v_key_num_cols"]:
            variant_key_txt += " \\\ " + var_txt_str
            col_idx = cols_for_str
        else:
            if col_idx > 0:
                variant_key_txt += " & "
            variant_key_txt += var_txt_str
            col_idx += cols_for_str
    
    variant_key_txt += '\end{tabular}'
    variant_key.text(setting_dict["v_key_x"], setting_dict["v_key_y"], variant_key_txt, verticalalignment='top', ha='center',
                     size=setting_dict["v_key_fontsize"], transform=variant_key.transAxes)

    return variant_key


def mark_up_special_chars(some_text):
    
    '''Mark up special characters for latex text.
    
    Args:
        | some_text (str):
    
    Returns:
        some_text: str
    '''
    
    some_text = some_text.replace("_", "\_")
    some_text = some_text.replace(">", "$>$")
    
    return some_text
        