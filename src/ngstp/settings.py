import pandas as pd

def get_setting_dict():

    '''
    Create a dictionary of settings to use in making the png. If the user wishes to change any of these default settings,
    then they should call this function to return the dictionary and change the values for the relevant keys.
    
    Returns
    -------
    setting_dict: dict
        settings for making the png.
    '''

    setting_dict = {}
    
    setting_dict["fig_width_inches"] = 11.69
    setting_dict["fig_height_inches"] = 8.27
    setting_dict["fig_dpi"] = 600
    setting_dict["fig_bottom"] = 0.1
    setting_dict["fig_top"] = 0.9
    setting_dict["fig_left"] = 0.125
    setting_dict["fig_num_rows"] = 17
    
    setting_dict["t_track_rows"] = 1
    setting_dict["c_track_rows"] = 3
    setting_dict["c_track_gap_rows"] = 2
    setting_dict["v_anns_top_rows"] = 2
    setting_dict["v_track_rows"] = 1
    setting_dict["v_anns_bot_rows"] = 2 
    setting_dict["v_key_rows"] = 2
    setting_dict["v_track_gap_rows"] = 1
    setting_dict["pd_track_rows"] = 4
    setting_dict["pd_track_gap_rows"] = 1
    
    setting_dict["title_1_fig_x"] = 0.5
    setting_dict["title_1_fig_y"] = 0.95
    setting_dict["title_2_ax_y"] = 0.5
    setting_dict["title_fontsize"] = 14
    
    setting_dict["c_track_y_axis_label"] = r'\noindent \textbf{Average}\\ \textbf{coverage}'
    setting_dict["c_track_fontsize"] = 10
    
    setting_dict["v_track_y_axis_label"] = r'\textbf{Exons \& variants}'
    setting_dict["v_track_vars_text_top"] = "Splice acceptor/donor \n (SA/D), initiator codon \n (IC), stop gained (SG),\n \& frameshift (F)"
    setting_dict["v_track_vars_text_bot"] = "Missense (M) \&\n inframe deletion (ID)"
    setting_dict["v_track_vars_text_top_x"] = -0.075 
    setting_dict["v_track_vars_text_top_y"] = 1.75
    setting_dict["v_track_vars_text_bot_x"] = -0.075
    setting_dict["v_track_vars_text_bot_y"] = -1.25
    setting_dict["v_track_arrow_style"] = "->"
    setting_dict["v_track_t_start"] = 1.7
    setting_dict["v_track_t_arrow_head"] = 0.95 
    setting_dict["v_track_t_inc"] = 0.4
    setting_dict["v_track_b_start"] = -0.9
    setting_dict["v_track_b_arrow_head"] = 0.075 
    setting_dict["v_track_b_inc"] = 0.4
    setting_dict["v_track_fontsize"] = 10
    setting_dict["v_track_var_abbrevs"] = {"missense_variant":"M","frameshift_variant":"F","stop_gained":"SG","splice_acceptor_variant":"SA",
                                       "splice_donor_variant":"SD","inframe_deletion":"ID","initiator_codon_variant":"IC"}
    setting_dict["v_track_vars_t_or_b"] = {"missense_variant":"B","frameshift_variant":"T","stop_gained":"T","splice_acceptor_variant":"T",
                                       "splice_donor_variant":"T","inframe_deletion":"B","initiator_codon_variant":"T"}
    setting_dict["v_track_merge_pixel_thresh"] = 3.5
    setting_dict["v_track_num_arrow_heights"] = 4
    
    setting_dict["v_key_num_cols"] = 4
    setting_dict["v_key_fontsize"] = 10
    setting_dict["v_key_x"] = 0.5 
    setting_dict["v_key_y"] = 1.175
    setting_dict["v_key_max_chars_per_col"] = 44 

    setting_dict["pd_track_y_label"] = r'\noindent\textbf{Protein \\ domains}'
    setting_dict["pd_track_y_label_fontsize"] = 10
    setting_dict["pd_track_y_label_x"] = -100.0
    setting_dict["pd_track_y_label_y"] = 0.5
    setting_dict["pd_track_y_label_va"] = 'center'
    setting_dict["pd_track_stripe_min_bases"] = 30
    setting_dict["pd_track_legend_fontsize"] = 10
    setting_dict["pd_track_legend_max_chars_per_row"] = 145
    setting_dict["pd_track_legend_bbox"] = (-0.025,-0.5,1.05,0.1)

    return setting_dict


def display_setting_dict(setting_dict):
    '''Pretty print the setting_dict.
    
    Parameters
    ----------
    setting_dict: dictionary
        settings for making the png.
    '''
    
    setting_dict_key_l = sorted(list(setting_dict.keys()))
    for key in setting_dict_key_l:
        print("{0}: {1}".format(key,setting_dict[key]))
    print("\n")
