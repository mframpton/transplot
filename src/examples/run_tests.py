import ngstp.ngs_transcript_plotter as ngstp
import ngstp.settings as s
import sys
import os


input_dir = os.path.abspath(os.path.join(os.getcwd(),"..","..","data","examples","input"))
#print(input_dir)
output_dir = os.path.abspath(os.path.join(os.getcwd(),"..","..","data","examples","plots"))

help(ngstp)
#sys.exit()

ngstp.make_protein_domain_color_file(
        protein_domain_file=os.path.join(input_dir, "APC_exoplot_domains_wt_overlaps.txt"),
        transcript_l= ["ENST00000457016"],
        database="Pfam", 
        sortby_col_l=["Start"],
        out_path=os.path.join(input_dir,"protein_domain_color.csv"))

ngstp.make_exon_coord_file(os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
                           "ENST00000457016",
                           os.path.join(input_dir , "APC_exon_coord.csv"))

ngstp.make_exon_coord_file(os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
                           "ENST00000457016",
                           os.path.join(input_dir , "APC_exon_coord_reverse.csv"))

setting_dict = s.get_setting_dict()
s.display_setting_dict(setting_dict)

ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["111"],
               [["543_A10"]],
               [os.path.join(input_dir, "APC_utrs.txt")],
               [os.path.join(input_dir, "APC_exon_coord.csv")],
               [os.path.join(input_dir, "APC_ENST00000457016_small.csv")],
               [os.path.join(input_dir, "APC_variants_CASES.txt")],
               [os.path.join(input_dir, "APC_exoplot_domains_wt_overlaps.txt")],
               os.path.join(input_dir, "protein_domain_color.csv"),
               setting_dict,
               os.path.join(output_dir , "APC_cases_make_png_1.png"))
#sys.exit()

'''
ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["110"],
               [["543_A10"]],
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "protein_domain_color.csv"),
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.2.png"))

ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["101"],
               [["543_A10"]],
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path(input_dir , "protein_domain_color.csv"),
               setting_dict,
               input_dir + "\\APC_cases.make_png.3.png")

ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["011"],
               [["543_A10"]],
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "protein_domain_color.csv"),
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.4.png"))

ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["010"],
               [["543_A10"]],
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "protein_domain_color.csv"),
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.5.png"))

ngstp.make_png(["ENST00000457016"],
               [r'\textbf{\textit{APC}}'],
               ["001"],
               [["543_A10"]],
               os.path.join(input_dir + "\\APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir + "\\APC_exon_coord_reverse.csv"),
               os.path.join(input_dir + "\\APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir + "\\APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir + "\\APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join( + "\\protein_domain_color.csv"),
               setting_dict,
               os.path.join( + "\\APC_cases.make_png.6.png"))


#ngstp.make_png(["ENST00000457016","ENST00000457016","ENST00000457016"],
#               [r'\textbf{\textit{APC1}}',r'\textbf{\textit{APC2}}',r'\textbf{\textit{APC3}}'],
#               ["100","100","100"],
#               [["543_A10"],["543_A10"],["543_A10"]],
#               os.path.join(input_dir + "\\APC_utrs.txt", input_dir + "\\APC_utrs_manual_reverse.txt", input_dir + "\\APC_utrs.txt"],
#            os.path.join(input_dir + "\\APC_exon_coord.csv", input_dir + "\\APC_exon_coord_reverse.csv", input_dir + "\\APC_exon_coord.csv"],
#               os.path.join(input_dir + "\\APC_ENST00000457016_small.csv", input_dir + "\\APC_ENST00000457016_small_reverse.csv", input_dir + "\\APC_ENST00000457016_small.csv"],
#               os.path.join(input_dir + "\\APC_variants_CASES.txt", input_dir + "\\APC_variants_CASES_reverse.txt", input_dir + "\\APC_variants_CASES.txt"],
#               os.path.join(input_dir + "\\APC_exoplot_domains_wt_overlaps.txt", input_dir + "\\APC_exoplot_domains_wt_overlaps.txt", input_dir + "\\APC_exoplot_domains_wt_overlaps.txt"],
#               input_dir + "\\protein_domain_color.csv",
#               setting_dict,
#               os.path.join(input_dir + "\\APC_cases.make_png.7.png"))


ngstp.make_png(["ENST00000457016","ENST00000457016"],
               [r'\textbf{\textit{APC1}}',r'\textbf{\textit{APC2}}'],
               ["010","010"],
               [["543_A10"],["543_A10"]],
               [input_dir + "\\APC_utrs.txt", input_dir + "\\APC_utrs_manual_reverse.txt"],
               [input_dir + "\\APC_exon_coord.csv", input_dir + "\\APC_exon_coord_reverse.csv"],
               [input_dir + "\\APC_ENST00000457016_small.csv", input_dir + "\\APC_ENST00000457016_small_reverse.csv"],
               [input_dir + "\\APC_variants_CASES.txt", input_dir + "\\APC_variants_CASES_reverse.txt"],
               [input_dir + "\\APC_exoplot_domains_wt_overlaps.txt", input_dir + "\\APC_exoplot_domains_wt_overlaps.txt"],
               input_dir + "\\protein_domain_color.csv",
               setting_dict,
               os.path.join(input_dir,"\\APC_cases.make_png.8.png"))


ngstp.make_png(["ENST00000457016","ENST00000457016","ENST00000457016"],
               [r'\textbf{\textit{APC1}}',r'\textbf{\textit{APC2}}',r'\textbf{\textit{APC3}}'],
               ["001","001","001"],
               [["543_A10"],["543_A10"],["543_A10"]],
               os.path.join(input_dir , "APC_utrs.txt"),
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt"),
               os.path.join(input_dir , "APC_utrs.txt"),
               os.path.join(input_dir , "APC_exon_coord.csv",
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_exon_coord.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
               os.path.join(input_dir , "APC_variants_CASES.txt"),
               input_dir + "\\APC_variants_CASES_reverse.txt", input_dir + "\\APC_variants_CASES.txt"],
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               input_dir + "\\APC_exoplot_domains_wt_overlaps.txt", input_dir + "\\APC_exoplot_domains_wt_overlaps.txt"],
               os.path.join(input_dir , "protein_domain_color.csv",
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.9.png")


ngstp.make_png(["ENST00000457016","ENST00000457016"],
               [r'\textbf{\textit{APC1}}',r'\textbf{\textit{APC2}}'],
               ["010","010"],
               [["543_A10"],["543_A10"]],
               os.path.join(input_dir , "APC_utrs.txt"),
               os.path.join(input_dir , "APC_exon_coord.csv"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_variants_CASES.txt"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "protein_domain_color.csv"),
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.10.png"))


ngstp.make_png(["ENST00000457016","ENST00000457016","ENST00000457016"],
               [r'\textbf{\textit{APC1}}',r'\textbf{\textit{APC2}}',r'\textbf{\textit{APC3}}'],
               ["010","010","010"],
               [["543_A10"],["543_A10"],["543_A10"]],
               os.path.join(input_dir , "APC_utrs.txt"),
               os.path.join(input_dir , "APC_utrs_manual_reverse.txt",
               os.path.join( input_dir , "APC_utrs.txt"),
               os.path.join(input_dir , "APC_exon_coord.csv"),
               os.path.join(input_dir , "APC_exon_coord_reverse.csv"),
               os.path.join( input_dir , "APC_exon_coord.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
               os.path.join(input_dir, "APC_ENST00000457016_small_reverse.csv"),
               os.path.join(input_dir , "APC_ENST00000457016_small.csv"),
               os.path.join(input_dir , "APC_variants_CASES.txt"),
               os.path.join(input_dir , "APC_variants_CASES_reverse.txt"),
               os.path.join(input_dir , "APC_variants_CASES.txt"),
               os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
              os.path.join( input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
              os.path.join(input_dir , "APC_exoplot_domains_wt_overlaps.txt"),
               os.path.join(input_dir , "protein_domain_color.csv"),
               setting_dict,
               os.path.join(input_dir , "APC_cases.make_png.11.png"))
'''