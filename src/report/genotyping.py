import os
import pandas as pd
import numpy as np
from sklearn.mixture import BayesianGaussianMixture
import warnings
import os
from matplotlib import pyplot as plt
import seaborn as sns
import src.templates as tmpl
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from matplotlib.collections import PolyCollection
from src.dtw_automata.overview import load_overview

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def calc_mae(true_alleles, preds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mae = np.mean([np.min([abs(i-lg) for lg in true_alleles]) for i in preds])
    return mae

def decode_alleles_complex(gmm_out_dict, df):
    if gmm_out_dict["is_hetero"]:
        df1 = df.loc[gmm_out_dict["group1"]]
        df2 = df.loc[gmm_out_dict["group2"]]
        mediangroup1 = []
        for col in df1.columns:
            if col!="reverse":
                mediangroup1.append(find_nearest(df1[col],np.median(df1[col])))
        mediangroup1_cnt = len(gmm_out_dict["group1"])
        mediangroup2 = []
        for col in df2.columns:
            if col!="reverse":
                mediangroup2.append(find_nearest(df2[col],np.median(df2[col])))
        mediangroup2_cnt = len(gmm_out_dict["group2"])
    else:
        df1 = df.loc[gmm_out_dict["group1"]]
        df2 = None
        mediangroup1 = []
        for col in df1.columns:
            if col!="reverse":
                mediangroup1.append(find_nearest(df1[col],np.median(df1[col])))
        mediangroup1_cnt = len(gmm_out_dict["group1"])
        mediangroup2 = '-'
        mediangroup2_cnt = '-'
    return (mediangroup1, mediangroup1_cnt, mediangroup2, mediangroup2_cnt)

def decode_alleles(gmm_out_dict):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if gmm_out_dict["is_hetero"]:
            mediangroup1 = find_nearest(gmm_out_dict["group1"],np.median(gmm_out_dict["group1"]))
            mediangroup1cnt = len(gmm_out_dict["group1"])
            mediangroup2 = find_nearest(gmm_out_dict["group2"],np.median(gmm_out_dict["group2"]))
            mediangroup2cnt = len(gmm_out_dict["group2"])
        else:
            mediangroup1 = find_nearest(gmm_out_dict["group1"],np.median(gmm_out_dict["group1"]))
            mediangroup1cnt = len(gmm_out_dict["group1"])
            mediangroup2 = '-'
            mediangroup2cnt = '-'
    return (mediangroup1,mediangroup1cnt,mediangroup2,mediangroup2cnt)

def filter_out_complex(df,cols,STD_COEFF):
    if len(df)<=5:
        return df
    
    idxes = []
    for col in cols:
        mean = np.mean(df[col])
        std = np.std(df[col])
        idxes.append([idx for idx,i in enumerate(df[col]) if i < (mean-STD_COEFF*std) or i > (mean+STD_COEFF*std)])
    
    idxes = np.unique([item for sublist in idxes for item in sublist])
    
    df.drop(idxes, inplace=True)
    
    return df

def filter_out(values,STD_COEFF):
    if len(values)<=5:
        return values

    mean = np.mean(values)
    std = np.std(values)
    out =  [i for i in values if i >= (mean-STD_COEFF*std) and i <= (mean+STD_COEFF*std)]
    return out

def run_genotyping_overview(overview, locus_path, config, muscle_path):
    if overview is None:
        overview_path, overview = load_overview(locus_path)

    unfilt_values, unfilt_basecalls = [], []
    if 'r_seq_start' in overview.columns and 'l_seq_end' in overview.columns:
        use_basecalls = True
    else:
        use_basecalls = False
    for row in overview.itertuples():
        if row.saved:
            unfilt_values.append(row.results)
            if use_basecalls:
                unfilt_basecalls.append(row.r_seq_start-row.l_seq_end)
    
    values = filter_out(unfilt_values, config["std_filter"])     
    gmm_out_dict = run_genotyping(values,config["min_weight"])

    if use_basecalls:
        basecalls = filter_out(unfilt_basecalls, config["std_filter"]) 
        gmm_out_dict_bc = run_genotyping(basecalls,config["min_weight"])  

    alleles = decode_alleles(gmm_out_dict)
    if use_basecalls:
        alleles_bc = decode_alleles(gmm_out_dict_bc) 

    final_preds_file = os.path.join(locus_path,"predictions","alleles.csv")
    with open(final_preds_file,"w") as f: 
        f.write("WarpSTR_allele1,WarpSTR_allele1_freq,WarpSTR_allele2,WarpSTR_allele2_freq,"\
                                                      "basecall_allele1,basecall_allele1_freq,"\
                    "basecall_allele2,basecall_allele2_freq\n")
        f.write("{a1},{a1f},{a2},{a2f},".format(a1=alleles[0],a1f=alleles[1],a2=alleles[2],a2f=alleles[3]))
        if use_basecalls:
            f.write("{b1},{b1f},{b2},{b2f}".format(b1=alleles_bc[0],b1f=alleles_bc[1],\
                                                         b2=alleles_bc[2],b2f=alleles_bc[3]))
    print(f"Allele lengths as given by WarpSTR: {alleles[0]},{alleles[2]}, frequency: {alleles[1]},{alleles[3]}")
    if use_basecalls:
        print(f"Allele lengths as given by basecall: {alleles_bc[0]},{alleles_bc[2]}, frequency: {alleles_bc[1]},{alleles_bc[3]}")
    
    if config['visualize']:
        img_path = os.path.join(locus_path,tmpl.SUMMARY_SUBDIR,"alleles.svg")
        vals = (gmm_out_dict["group1"],gmm_out_dict["group2"])
        if use_basecalls:
            bvals = (gmm_out_dict_bc["group1"],gmm_out_dict_bc["group2"])
            plot_clustering_preds(img_path, vals, bvals, alleles, alleles_bc)
        else:
            plot_clustering_preds(img_path, vals, None, alleles, None)

    if config['msa']:
        path_our = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,"sequences")
        msa_out, trcalls = call_muscle(muscle_path, path_our)

        if use_basecalls:
            path_bc = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,"basecalls")
            msa_out_bc, basecalls = call_muscle(muscle_path, path_bc)
        if gmm_out_dict["is_hetero"]:
            msa_out_group1, msa_out_group2 = groups_from_pred(muscle_path, path_our, gmm_out_dict, trcalls)
        
        if gmm_out_dict_bc["is_hetero"] and use_basecalls:
            msa_out_group1_bc, msa_out_group2_bc = groups_from_pred(muscle_path, path_bc, gmm_out_dict_bc, basecalls)

def run_genotyping_complex(locus_path, config,locus, df):
    if df is None:
        inpath = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,tmpl.COMPLEX_SUBDIR,"complex_repeat_units.csv")
        df = pd.read_csv(inpath,index_col=0)
    
    cols = [col for col in df.columns if col != "reverse"]
    if len(cols)<2:
        return

    df = filter_out_complex(df, cols, config["std_filter"])
    X = np.array(df[cols]).reshape(-1, len(cols))
    model_bayes = BayesianGaussianMixture(weight_concentration_prior=0.25, covariance_type='tied', n_components=2,n_init=5,max_iter = 1000).fit(X)

    out = {}
    if any(x < config["min_weight"] for x in model_bayes.weights_):
        out["is_hetero"] = False
        out["group1"] = df.index
        out["group2"] = []
        out["predictions"] = None
    else:
        preds = model_bayes.predict(X)
        out["is_hetero"] = True
        out["group1"] = [idx for idx,g in zip(df.index,preds) if g==0]
        out["group2"] = [idx for idx,g in zip(df.index,preds) if g==1]
        out["predictions"] = preds    

    alleles = decode_alleles_complex(out,df)
    if out["predictions"] is not None:
        df["allele"] = preds

    img_path = os.path.join(locus_path,tmpl.SUMMARY_SUBDIR,"complex_genotypes.svg")
    fig, axes = plt.subplots(nrows=len(cols), ncols=1, figsize=(8, 6*len(cols)))
    for idx,col in enumerate(cols):
        ex_df = df
        val1 = alleles[0][idx]
        val2 = alleles[2][idx]
        name = "repeat numbers: "+str(val1)+","+str(val2)
        ex_df[name] = ""
        axes[idx] = sns.violinplot(data=ex_df, x=name,y=col,hue="allele",orient="vertical", split=False,scale="count",whis=np.inf,inner=None,ax=axes[idx])
        plt.setp(axes[idx].collections, alpha=.3)
        first = [r for r in axes[idx].get_children() if type(r)==PolyCollection]
        c1 = first[0].get_facecolor()[0]
        c2 = first[1].get_facecolor()[0]
        if val1!="-":
            axes[idx].axhline(y=val1, color=c1, linestyle='--')
        if val2!="-":
            axes[idx].axhline(y=val2, color=c2, linestyle='--')
        axes[idx] = sns.stripplot(data=ex_df, x=name,y=col,hue="allele",orient="vertical",dodge=True,size=6,alpha=0.8,jitter=0.3,ax=axes[idx])
        axes[idx].get_legend().remove()
    
    #fig.suptitle('Clustered predictions of repeat units \n Config sequence: '+locus["sequence"])
    plt.savefig(img_path,bbox_inches='tight',format="svg")
    plt.close()  

def run_genotyping(vals,MIN_WEIGHT):
    out = {}
    if len(np.unique(vals))==1:
        out["is_hetero"] = False
        out["group1"] = vals
        out["group2"] = []
    
    else:
        X = np.array(vals).reshape(-1, 1)
        model_bayes = BayesianGaussianMixture(weight_concentration_prior=0.25, covariance_type='tied', n_components=2,n_init=5,max_iter = 1000).fit(X)

        if any(x < MIN_WEIGHT for x in model_bayes.weights_):
            out["is_hetero"] = False
            out["group1"] = vals
            out["group2"] = []

        else:
            preds = model_bayes.predict(X)
            out["is_hetero"] = True
            out["group1"] = [i for i,g in zip(vals,preds) if g==0]
            out["group2"] = [i for i,g in zip(vals,preds) if g==1]
            out["predictions"] = preds    
    return out

def call_muscle(muscle_path, path):
    if os.path.exists(muscle_path) is False:
        raise ValueError(f'Could not load the MUSCLE tool from {muscle_path}')
    calls_file = os.path.join(path,"all.fasta")
    msa_out = os.path.join(path,"msa_all.fasta")
    calls = [(record.id, record.seq) for record in SeqIO.parse(calls_file, "fasta")]
    muscle_cline = MuscleCommandline(muscle_path, input=calls_file, out=msa_out,verbose=False,quiet=True)
    os.system(str(muscle_cline))
    return msa_out,calls

def groups_from_pred(muscle_path, path_out, gmm_out_dict, calls):
    groups = gmm_out_dict["predictions"]
    group1_seqs = [i for i,g in zip(calls,groups) if g==0]
    group2_seqs = [i for i,g in zip(calls,groups) if g==1]
            
    trcalls_file_group1 = os.path.join(path_out,"group1.fasta")
    store_fasta(trcalls_file_group1,group1_seqs)

    msa_out_group1 = os.path.join(path_out,"msa_group1.fasta")
    muscle_cline = MuscleCommandline(muscle_path, input=trcalls_file_group1, out=msa_out_group1,verbose=False,quiet=True)
    os.system(str(muscle_cline))

    trcalls_file_group2 = os.path.join(path_out,"group2.fasta")
    store_fasta(trcalls_file_group2,group2_seqs)
         
    msa_out_group2 = os.path.join(path_out,"msa_group2.fasta")
    muscle_cline = MuscleCommandline(muscle_path, input=trcalls_file_group2, out=msa_out_group2,verbose=False,quiet=True)
    os.system(str(muscle_cline))
    return msa_out_group1,msa_out_group2

def store_fasta(filename,fasta):
    with open (filename,"w") as file:
        for seqid,seqrec in fasta:
                file.write(">"+seqid+"\n")
                file.write(str(seqrec)+"\n\n")

def plot_clustering_preds(img_path, vals, bvals, alleles, alleles_bc):
    if bvals is not None:
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 7),sharey=True,gridspec_kw={'width_ratios': [1, 1]})
        axes[0] = sns.violinplot(data=vals, orient="vertical", split=True,scale="count",whis=np.inf,inner=None, color=".8",width=0.8,scale_hue=False,ax=axes[0])
        axes[0] = sns.stripplot(data=vals, orient="vertical",size=6,alpha=0.6,jitter=0.3,palette="bright",ax=axes[0])
        axes[0].title.set_text('WarpSTR')
        axes[0].set(xlabel='Alleles of length {p1},{p2}'.format(p1=alleles[0],p2=alleles[2]))
        axes[0].set(ylabel='Filtered predictions')
        axes[1] = sns.violinplot(data=bvals, orient="vertical", split=True,scale="count",whis=np.inf,inner=None, color=".8",width=0.8,scale_hue=False,ax=axes[1])
        axes[1] = sns.stripplot(data=bvals, orient="vertical",size=6,alpha=0.6,jitter=0.3,palette="bright",ax=axes[1])
        axes[1].title.set_text('Basecalled sequences')
        axes[1].set(xlabel='Alleles of length {p1},{p2}'.format(p1=alleles_bc[0],p2=alleles_bc[2]))
        axes[1].set(ylabel='Filtered predictions')
        plt.savefig(img_path,bbox_inches='tight',format="svg")
        plt.close()
        return
    
    fig = plt.figure(figsize=(16, 7))
    fig.suptitle('WarpSTR', fontsize=20)
    plt.xlabel('Alleles of length {p1},{p2}'.format(p1=alleles[0],p2=alleles[2]), fontsize=18)
    plt.ylabel('Filtered predictions', fontsize=16)
    sns.violinplot(data=vals, orient="vertical", split=True,scale="count",whis=np.inf,inner=None, color=".8",width=0.8,scale_hue=False)
    sns.stripplot(data=vals, orient="vertical",size=6,alpha=0.6,jitter=0.3,palette="bright")
    plt.savefig(img_path,bbox_inches='tight',format="svg")
    plt.close()
