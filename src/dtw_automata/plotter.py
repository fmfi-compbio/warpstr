import os
import bisect
import seaborn as sns
from matplotlib import pyplot as plt
import src.templates as tmpl


def plot_collapsed(df_overview, locus_path):
    """
    Plot predictions as per repeat unit
    """
    cols = [col for col in df_overview.columns if col != "reverse"]
    out_path = os.path.join(
        locus_path, tmpl.SUMMARY_SUBDIR, "collapsed_predictions.svg")
    fig = plt.figure(figsize=(12, 6))
    sns.violinplot(data=df_overview[cols], orient="vertical", split=True,
                   scale="count", inner=None, color=".8", width=0.8, scale_hue=False)
    sns.stripplot(data=df_overview[cols], orient="vertical",
                  size=6, alpha=0.6, jitter=0.3, palette="bright")
    plt.ylabel('Predicted repeat numbers')
    fig.suptitle('Predictions of repeat units')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()

    df_template = df_overview[df_overview["reverse"] == False]
    df_reverse = df_overview[df_overview["reverse"] == True]
    out_path = os.path.join(locus_path, tmpl.SUMMARY_SUBDIR,
                            "collapsed_predictions_strand.svg")
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(
        16, 7), sharey=True, gridspec_kw={'width_ratios': [1, 1]})
    axes[0] = sns.violinplot(data=df_template[cols], orient="vertical", split=True,
                             scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[0])
    axes[0] = sns.stripplot(data=df_template[cols], orient="vertical",
                            size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[0])
    axes[0].set(xlabel='Template strand')
    axes[0].set(ylabel='Predicted repeat numbers')
    axes[1] = sns.violinplot(data=df_reverse[cols], orient="vertical", split=True,
                             scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[1])
    axes[1] = sns.stripplot(data=df_reverse[cols], orient="vertical",
                            size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[1])
    axes[1].set(xlabel='Reverse strand')
    fig.suptitle('Predictions of repeat units split by strand')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()


def save_warp_img(out_warp_path, flank_length, warper, resc_warper, reverse, read_name):
    """
    Save warping and rescaling image
    """
    out_path = os.path.join(out_warp_path, read_name+"_" +
                            str(reverse)+"_"+str(len(resc_warper.seq))+".png")
    start, end, resc_start, resc_end = delineate_tr_region(
        flank_length, warper.trace, resc_warper.trace, (warper.offset, resc_warper.offset))
    note = None
    if (end-start) > 2000 or (resc_end-resc_start) > 2000:
        note = "*Truncated to 2000 as original length was " + \
            str(resc_end-resc_start)
        end = start + 2000
        resc_end = resc_start + 2000

    fig, axs = plt.subplots(2, figsize=(16, 4), sharey=True)
    axs[0].plot(warper.template[start:end], label='signal')
    axs[0].plot(warper.warped_signal[start:end], label='warped path')
    #axs[0].set_title(warper.seq)
    axs[0].set_ylabel('Norm. current')
    axs[0].legend()
    axs[1].plot(resc_warper.template[resc_start:resc_end],
                label='rescaled signal')
    axs[1].plot(resc_warper.warped_signal[resc_start:resc_end],
                label='warped path')
    #axs[1].set_title(resc_warper.seq)
    axs[1].set_ylabel('Norm. current')
    axs[1].legend()
    if note:
        plt.figtext(0.5, 0.01, note, ha="center", fontsize=12)
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()


def delineate_tr_region(flank_length, trace, resc_trace, offsets):
    """
    Delineates TR regions
    """
    boundary = flank_length - 20
    start = bisect.bisect_left(trace, boundary - offsets[0] - 5)  # -boundary
    maxstate = trace[-1]
    end = bisect.bisect_right(trace, maxstate - boundary)  # +boundary

    resc_start = bisect.bisect_left(
        resc_trace, boundary - offsets[1] - 5)  # -boundary
    resc_maxstate = resc_trace[-1]
    resc_end = bisect.bisect_right(
        resc_trace, resc_maxstate - boundary)  # +boundary
    return start, end, resc_start, resc_end


def plot_summaries(locus_path, df_out, config):
    """
    Plots summaries for results
    """
    df_out["dtw_cost2"] = df_out["dtw_cost2"].astype(float)
    df_out["results"] = df_out["results"].astype(float)
    df_out["reverse"] = df_out["reverse"].astype(float)

    df_saved = df_out[df_out['saved'] == 1]

    if 'r_seq_start' in df_out and 'l_seq_end' in df_out:
        df_out["bc_length"] = df_out["r_seq_start"]-df_out["l_seq_end"]
        use_basecall = True
    else:
        use_basecall = False
        df_out["bc_length"] = 0

    df_rev = df_out[(df_out['reverse'] == 1) & (df_out['saved'] == 1)]
    df_temp = df_out[(df_out['reverse'] == 0) & (df_out['saved'] == 1)]

    if config['visualize_strand']:
        out_path = os.path.join(
            locus_path, tmpl.SUMMARY_SUBDIR, "predictions_strand.svg")
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(
            16, 7), sharey=True, gridspec_kw={'width_ratios': [1, 1]})
        axes[0] = sns.violinplot(data=df_temp["results"], orient="vertical", split=True,
                                 scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[0])
        axes[0] = sns.stripplot(data=df_temp["results"], orient="vertical",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[0])
        axes[0].set(ylabel='Predicted allele lengths')
        axes[0].set(xlabel='Template strand')
        axes[1] = sns.violinplot(data=df_rev["results"], orient="vertical", split=True,
                                 scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[1])
        axes[1] = sns.stripplot(data=df_rev["results"], orient="vertical",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[1])
        axes[1].set(xlabel='Reverse strand')
        fig.suptitle('Predictions split by strand')
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()

        df_out.rename(columns={"A": "a", "B": "c"})
        df_rev.rename(columns={"results": "WarpSTR",
                      "bc_length": "Basecalling"}, inplace=True)
        df_temp.rename(columns={"results": "WarpSTR",
                       "bc_length": "Basecalling"}, inplace=True)
        df_rev_melt = df_rev.reset_index().melt(
            id_vars='read_name', value_vars=['WarpSTR', 'Basecalling'])
        df_temp_melt = df_temp.reset_index().melt(
            id_vars='read_name', value_vars=['WarpSTR', 'Basecalling'])

        out_path = os.path.join(
            locus_path, tmpl.SUMMARY_SUBDIR, "predictions_strand_with_basecalls.svg")
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(
            16, 7), sharey=True, gridspec_kw={'width_ratios': [1, 1]})
        axes[0] = sns.violinplot(data=df_temp_melt, x="variable", y="value", scale="count",
                                 inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[0])
        axes[0] = sns.stripplot(data=df_temp_melt, x="variable", y="value",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[0])
        axes[0].set(ylabel='Predicted allele lengths')
        axes[0].set(xlabel='Template strand')
        axes[1] = sns.violinplot(data=df_rev_melt, x="variable", y="value", scale="count",
                                 inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[1])
        axes[1] = sns.stripplot(data=df_rev_melt, x="variable", y="value",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[1])
        axes[1].set(xlabel='Reverse strand')
        axes[1].set(ylabel='Predicted allele lengths')
        fig.suptitle('Predictions split by strand')
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()

    if config['visualize_phase']:
        out_path = os.path.join(
            locus_path, tmpl.SUMMARY_SUBDIR, "predictions_phase.svg")
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(
            16, 7), sharey=True, gridspec_kw={'width_ratios': [1, 1]})
        axes[0] = sns.violinplot(data=df_saved["orig"], orient="vertical", split=True,
                                 scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[0])
        axes[0] = sns.stripplot(data=df_saved["orig"], orient="vertical",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[0])
        axes[0].set(ylabel='Predicted allele lengths')
        axes[0].set(xlabel='First phase')
        axes[1] = sns.violinplot(data=df_saved["results"], orient="vertical", split=True,
                                 scale="count", inner=None, color=".8", width=0.8, scale_hue=False, ax=axes[1])
        axes[1] = sns.stripplot(data=df_saved["results"], orient="vertical",
                                size=6, alpha=0.6, jitter=0.3, palette="bright", ax=axes[1])
        axes[1].set(xlabel='Second phase')
        fig.suptitle('Predictions split by phase')
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()

    if config['visualize_cost']:
        out_path = os.path.join(
            locus_path, tmpl.SUMMARY_SUBDIR, "predictions_cost.svg")
        axis = sns.scatterplot(data=df_saved, x="results", y="dtw_cost2")
        axis.set(xlabel='Predicted allele lengths',
                 ylabel='Calculated state-wise cost')
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()
