import io
import pandas as pd
import os
import src.templates as tmpl
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline


def load_reference(locus_path):
    """
    Function loading template and reverse flanks
    :param locus_path: Path to the extracted data for selected locus
    :returns (lft,rft): left and right flanking template sequences
    :returns (lfr,rfr): left and right flanking reverse sequences
    """
    path = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_FLANKS)
    lft,rft,lfr, rfr = "","","",""
            
    if os.path.exists(path):   
        with open(path, "r") as f:
            colnames = f.readline().split(',')
            seqid = 0
            for idx,i in enumerate(colnames):
                if i.rstrip()=="sequence":
                    seqid = idx
            
            for i in range(4):
                f.readline()
                
            template = f.readline().split(',')[seqid].rstrip()
        
    else:
        print("Error: loading of refernce failed - "+path+" not Found")
    
    if len(template)<1:
        print("Error when loading reference sequence during reporting")
        return False
    
    return template

def load_signal(path):
    sig = []
    if os.path.exists(path):
        with open(path, "r") as f:
            sig = [float(i.strip()) for i in f.readlines()]
    else:
        print("Error: loading of exp signal failed - "+path+" not Found")
    return sig
        
def load_reference_signal(locus_path,flanksize=10):
    path_temp_left = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[0]+".txt")
    path_temp_right = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[1]+".txt")
    path_temp_rep = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[4]+".txt")
    
    temp_left = load_signal(path_temp_left)
    temp_right = load_signal(path_temp_right)
    temp_repeat = load_signal(path_temp_rep)
    
    if len(temp_left)>0 and len(temp_repeat)>0 and len(temp_right)>0:
        temp_signal = temp_left[-flanksize:]+temp_repeat+temp_right[:flanksize]
    else:
        print("Error when loading reference template signal during reporting")
        temp_signal = []
    
    path_rev_left = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[2]+".txt")
    path_rev_right = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[3]+".txt")
    path_rev_rep = os.path.join(locus_path,tmpl.LOCUS_INFO_SUBDIR,tmpl.LOCUS_NAMES[5]+".txt")
    
    rev_left = load_signal(path_rev_left)
    rev_right = load_signal(path_rev_right)
    rev_repeat = load_signal(path_rev_rep)
    
    if len(rev_left)>0 and len(rev_repeat)>0 and len(rev_right)>0:
        rev_signal = rev_right[-flanksize:]+rev_repeat+rev_left[:flanksize]
    else:
        print("Error when loading reference template signal during reporting")
        rev_signal = []
    
    return temp_signal,rev_signal

def get_summary_overview(df_overview):
    """
    :param: pandas dataframe of overview file
    returns summary
    """
    df_out = df_overview.groupby(['sample','run_id'])['sample'].count()
    samples_df = (pd.concat([df_out.to_frame(),df_out.sum(level=0)\
                           .to_frame().assign(claim_type= "TOTAL:")\
                           .set_index('claim_type', append=True)]).sort_index())
    df_saved = df_overview[df_overview['saved']==1]
    
    new = []
    total = 0
    for r,i in samples_df.index.values:
        if i!="TOTAL:":
            sz = len(df_saved[df_saved["run_id"]==i])
            new.append(sz)
            total += sz 
        else:
            new.append(total)
            total = 0
    samples_df.rename(columns = {'sample':'mapped reads'}, inplace = True) 
    samples_df['extracted']=new
    return samples_df

def load_textfile(path):
    s = ''
    with open(path,"r") as f: 
        for i in f.readlines():
            s += i.strip()
    return s

def get_scripts(path):
    script = '<script>'+load_textfile(path)+'</script>'
    return script

def get_style(path):
    style = '<style>'+load_textfile(path)+'</style>'
    return style

def load_overview(path):
    df = pd.read_csv(path)
    df.set_index('read_name',inplace=True)
    df.columns = df.columns.map(str)
    return df

def create_html_table(table_id,df,klass):
    str_io = io.StringIO()
    df.to_html(buf=str_io, classes=klass, table_id=table_id)
    html_table = str_io.getvalue()
    return html_table

def get_head(name):
    head = '<!DOCTYPE HTML><html><head>  <meta charset="UTF-8">    <meta name="description"'
    head += 'content="'+name+'" /><title>'+name+'</title>'
    return head


def get_results_overview(df_overview,locus_path,results):
    
    basic = []
    idx = 0
    for row in df_overview.itertuples():
        if row.saved:
            bq_len = row.r_seq_start-row.l_seq_end
            i,j = results[idx][0],results[idx][1]
            basic.append((row.Index,row.sample,row.reverse,bq_len,i,j,"test","test","test"))
            idx = idx + 1
            
    df = pd.DataFrame(basic,columns=("read_name","sample","reverse", "Basecalling", "DTW", "RescDTW","Basecalling_seq", "DTW_seq", "RescDTW_seq"))
    
    df.set_index('read_name',inplace=True)
    return df

def get_ref_svgs(locus_path,imgpath):
    temp_signal,rev_signal = load_reference_signal(locus_path)
        
    plt.figure(figsize=[7, 3])
    plt.plot(temp_signal,'-o')
    #plt.axis('off')
    plt.title("template ref signal")
    plt.gca().set_position([0, 0, 1, 1])
    plt.axvline(x = 9, color = 'r')
    plt.axvline(x = len(temp_signal)-10, color = 'r')
    plt.savefig(imgpath,bbox_inches='tight')
    plt.close()     
    svg1 = load_textfile(imgpath)  
    
    plt.figure(figsize=[7, 3])
    plt.plot(rev_signal,'-o')
    #plt.axis('off')
    plt.title("reverse ref signal")
    plt.gca().set_position([0, 0, 1, 1])
    plt.axvline(x = 9, color = 'r')
    plt.axvline(x = len(rev_signal)-10, color = 'r')
    plt.savefig(imgpath,bbox_inches='tight')
    plt.close()     
    svg2 = load_textfile(imgpath)
    os.remove(imgpath) 
    return svg1,svg2


def filter_out(values,STD_COEFF=2):
    if len(values)<=5:
        return values,range(len(values))

    mean = np.mean(values)
    std = np.std(values)
    out =  [i for i in values if i >= (mean-STD_COEFF*std) and i <= (mean+STD_COEFF*std)]
    idxes = [idx for idx,i in enumerate(values) if i >= (mean-STD_COEFF*std) and i <= (mean+STD_COEFF*std)]
    #print("removed",len(idxes)-len(values),"from",len(values))
    return out,idxes

def call_muscle(muscle_path, locus_path,flag):
    calls_file = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,flag+".fasta")
    msa_out = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,flag+"_msa.fasta")
    calls = [(record.id, record.seq) for record in SeqIO.parse(calls_file, "fasta")]
    muscle_cline = MuscleCommandline(muscle_path, input=calls_file, out=msa_out,verbose=False,quiet=True)
    os.system(str(muscle_cline))
    return msa_out,calls

def plot_clustering_preds(gmm_out_dict, gmm_out_dict_bc, allele_lengths, mae_basecalls, mae_our, img_path, vals, bvals, error):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 7),sharey=True,gridspec_kw={'width_ratios': [1, 1]})
    axes[0] = sns.violinplot(data=vals, orient="vertical", split=True,scale="count",whis=np.inf,inner=None, color=".8",width=0.8,scale_hue=False,ax=axes[0])
    axes[0] = sns.stripplot(data=vals, orient="vertical",size=6,alpha=0.6,jitter=0.3,palette="bright",ax=axes[0])
    for i in allele_lengths:
        axes[0].axhline(i, color ='green', lw = 2, alpha = 0.75,linestyle ="--")
    axes[0].title.set_text('WarpSTR')
    axes[0].set(ylabel='Lengths for alleles '+prepare_plot_label(gmm_out_dict))
    axes[0].set(xlabel='MAE:'+str(error[0])+"("+str(mae_our)+")")
    axes[1] = sns.violinplot(data=bvals, orient="vertical", split=True,scale="count",whis=np.inf,inner=None, color=".8",width=0.8,scale_hue=False,ax=axes[1])
    axes[1] = sns.stripplot(data=bvals, orient="vertical",size=6,alpha=0.6,jitter=0.3,palette="bright",ax=axes[1])
    for i in allele_lengths:
        axes[1].axhline(i, color ='green', lw = 2, alpha = 0.75,linestyle ="--")
    axes[1].title.set_text('Basecalled sequences')
    axes[1].set(ylabel='Lengths for alleles '+prepare_plot_label(gmm_out_dict_bc))   
    axes[1].set(xlabel='MAE:'+str(error[1])+"("+str(mae_basecalls)+")")
    plt.savefig(img_path,bbox_inches='tight',format="svg")
    plt.close()     
    svg = load_textfile(img_path)  
    svg_string = svg + "<br>"
    return svg_string


def get_index_html(all_loci):
    html = '<a name="Index"><div id="Index">'
    html += "<h1>Index</h1><ul>"
    for locus,_ in all_loci:
        html += '<li><a href="#'+locus['name']+'">'+locus['name']+'</a></li>'
    html += '</ul></div></a>'
    return html

def get_locus_header(locus):
    html = "<hr>"
    html += '<a name="'+locus['name']+'"><div id="'+locus['name']+'">'
    html += "<h1>"+locus['name']+"</h1>"
    html += "</div></a>"
    html += "<h4>Ref locus: "+locus['noting']+"</h4>"
    html += "<h4>Coding sequence: "+locus['sequence']+"</h4>"
    return html

def parse_vcf_string(locus,locus_path):
    svg_string = ""
    allele_lengths = []
    if "vcf_output" in locus:
        if type(locus["vcf_output"]) is int:
            allele_lengths.append(locus["vcf_output"])
        else:
            allele_lengths = [int(i) for i in locus["vcf_output"].split(",")]
        svg_string += "<h4>*Showing other method prediction as the green dashed line</h4><br>"
    else:
        allele_lengths.append(len(load_reference(locus_path)))
        svg_string += "<h4>*No VCF output - Showing reference allele length as the green dashed line</h4><br>"
    return allele_lengths,svg_string

def get_msa_html(msa_file,name,title,allele):
    with open(msa_file,"r") as msaf:
        msa = msaf.read()
    msa_div_id = "msa_"+name
    fasta_id = "fasta-file-"+name
    s = '<div id='+msa_div_id+'>Loading Multiple Alignment...</div>'        
    s += '<pre style="display: none" id='+fasta_id+'>'+msa+'</pre>'
    msa_js = 'var fasta = document.getElementById("'+fasta_id+'").innerText; var m = msa({el: document.getElementById("'+msa_div_id+'"),seqs: msa.io.fasta.parse(fasta), colorscheme: {"scheme": "nucleotide"},bootstrapMenu: true,menu: "small",vis: {conserv: false,metaIdentity: true,seqlogo: true}});m.render();' 
    
    html_msa = wrap_collapsible("Show alignment visualization for "+title+" output"+allele,s,active=False) + "<br>"
    script_msa = '<script>'+msa_js+'</script>'
    return html_msa+script_msa

def prepare_plot_label(gmm_out_dict):
    a = ""
    if len(gmm_out_dict["group1"])>0:
        a+=str(np.median(gmm_out_dict["group1"]))
    if len(gmm_out_dict["group2"])>0:
        a+=","+str(np.median(gmm_out_dict["group2"]))
    return a

def wrap_collapsible(title,data,active=False):
    if active:
        button_tag_start1 = '<button type="button" class="collapsible active">'
    else:
        button_tag_start1 = '<button type="button" class="collapsible">'
    
    if active:
        button_tag_start2 = '</button><div class="content" style="display: block;">'
    else:
        button_tag_start2 = '</button><div class="content">'
    button_end='</div>'
    s = button_tag_start1+title+button_tag_start2+data+button_end
    return s

def create_report(main_out_path, tr_results_seq,locus,locus_path,src_path):
    select_box_path = os.path.join(src_path,"selectbox.html")
    select_box = load_textfile(select_box_path)

    reports_path = os.path.join(main_out_path,tmpl.REPORTS_SUBDIR)
    output_file = os.path.join(reports_path,tmpl.REPORT_FILENAME.format(locus=locus['name']))
    
    reference = load_reference(locus_path)
    overview_file = os.path.join(locus_path,tmpl.OVERVIEW_NAME)
    
    df_overview = load_overview(overview_file)
    summary_df = get_summary_overview(df_overview)
    
    filtering = False
    meanres = False

    #results_df = filter_out(df_overview, locus_path, tr_results)
    #print("filtering",results_df)
    #if filtering:
    #     results_df = filter_out(df_overview, locus_path, tr_results)
    #elif meanres:
    #    results_df = meanresults(df_overview, locus_path, tr_results)
    #    print("meanres",results_df)
    #else:
    tr_results = []
    for seq,seq_resc in tr_results_seq:
        tr_results.append((len(seq),len(seq_resc)))

    results_df = get_results_overview(df_overview, locus_path, tr_results)
    
    try:
        with open(output_file,"w+") as f:
            f.write(get_head(tmpl.REPORT_FILENAME.format(locus=locus['name'])))
            f.write("<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>")
            f.write(get_style(os.path.join(src_path,"style.css")))
            f.write('</head><body>')

            f.write("<h2>Report for "+locus['name']+"</h2>")
            f.write(create_html_table('summaryTable',summary_df,'summarytable'))      
            f.write("<h2>Locus reference: "+locus['noting']+"</h2>")
            f.write("<h2>"+reference+"</h2>")
            f.write("<h2>TR length: "+str(len(reference))+"</h2>")
            f.write("<h2 id='input_sequence'>Our input sequence: "+locus['sequence']+"</h2>")
            f.write('<br><br>')
            f.write(select_box)
            #f.write("<p id='textareabasic1'>Click first to show sequences</p>")
            #f.write("<p id='textareabasic2'>Click first to show sequences</p>")
            #f.write("<p id='textareabasic3'>Click first to show sequences</p>")
            f.write("<div id='boxplot_div'></div>")
            f.write('<br><br>')
            f.write(create_html_table('myTable',results_df,'myTable'))
            #f.write(create_html_table('othertable',gmm_results,'myTable'))
            
            #tr_calls_reports = os.path.join(locus_path,tmpl.PREDICTIONS_SUBDIR,tmpl.REPORTS_SUBDIR)
            #automat1_path = os.path.join(tr_calls_reports,"automat_simple_template.svg")
            #print(automat1_path,flush=True)
            #if os.path.exists(automat1_path):
            #    f.write('<br><br>')
            #    f.write("<div id='simple_states'>")
            #    f.write(load_textfile(automat1_path))
            #    f.write("</div>")
            #
            #automat2_path = os.path.join(tr_calls_reports,"automat_kmer_template.svg")
            #if os.path.exists(automat2_path):
            #    f.write('<br><br>')
            #    f.write("<div id='simple_states'>")
            #    f.write(load_textfile(automat2_path))
            #    f.write("</div>")
            
            f.write(get_scripts(os.path.join(src_path,"filter_table.js")))
            #f.write(get_scripts(os.path.join(src_path,"prepare_selectbox.js")))
            f.write(get_scripts(os.path.join(src_path,"boxplots.js")))
            f.write("</body></html>")    
            print(tmpl.REPORT_SUCCESS.format(path=output_file))    
    except:
        print(tmpl.REPORT_FAILURE.format(path=output_file))  
