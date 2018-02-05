import pandas as pd

'''
Functions for the UTRs file.
'''

def get_utr_df(utr_file, strand, transcript):
    
    '''
    Read the UTR information into a DataFrame.
    
    Parameters
    ----------
    utr_file: str
        path to file containing the utr information.
    strand: str
        whether the transcript is on the positive or negative strand.
    transcript: str
        Ensemble transcript ID.
    
    Returns
    -------
    utr_df: DataFrame
    '''
    
    utr_df = pd.DataFrame([])
    utr_df_chunker = pd.read_csv(utr_file,sep=",",chunksize=1000)
    for utr_df_chunk in utr_df_chunker:
        utr_df_chunk = utr_df_chunk[utr_df_chunk["Ensembl Transcript ID"] == transcript]
        utr_df = pd.concat([utr_df, utr_df_chunk], ignore_index=True)
        
    utr_df = utr_df[~utr_df["5' UTR Start"].isnull() | ~utr_df["3' UTR Start"].isnull()]
    #Convert the 0-based start coordinates to 1-based.
    utr_df["5' UTR Start"] = utr_df["5' UTR Start"] + 1
    utr_df["3' UTR Start"] = utr_df["3' UTR Start"] + 1
    
    #Get utr_df into the desired form.
    utr_df["utr"] = utr_df["5' UTR Start"].isnull()
    utr_df["utr"].replace(to_replace={False:"5'", True:"3'"}, inplace=True)
    utr_df["5' UTR Start"].fillna(utr_df["3' UTR Start"], inplace=True)
    utr_df["5' UTR End"].fillna(utr_df["3' UTR End"], inplace=True)
    utr_df = utr_df[["utr", "5' UTR Start", "5' UTR End"]]
    utr_df.rename(columns={"5' UTR Start":"start_bp", "5' UTR End":"end_bp"}, inplace=True)

    if strand == "-":
        utr_df.rename(columns={"start_bp":"end_bp", "end_bp":"start_bp", "start_tp":"end_tp", "end_tp":"start_tp"}, inplace=True)

    return utr_df
