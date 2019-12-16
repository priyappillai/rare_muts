import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pickle
from scipy.stats import beta, kstest
from scipy.special import btdtr


pops = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas']
namepops = ["African", "Latino", "Ashkenazi Jewish", "East Asian", "European (Finnish)", "European (non-Finnish)", "South Asian"]
acpops = ['AC_%s=' % pop for pop in pops]
anpops = ['AN_%s=' % pop for pop in pops]
afpops = ['AF_%s=' % pop for pop in pops]

variant_types = ["splice", "start_lost", "stop_lost", "stop_gained", "frameshift", 
                 "inframe_insertion", "inframe_deletion", "missense_variant", "UTR", 
                 "synonymous", "stop_retained", "intron", "downstream", "upstream", 
                 "protein_altering", "non_coding", "intergenic", "regulatory_region", 
                 "incomplete_terminal_codon", "mature_miRNA", "protein_coding", 
                 "coding_sequence"]

nonsyn_coding_variants = ["start_lost", "stop_lost", "stop_gained", "frameshift", 
                          "inframe_insertion", "inframe_deletion", "missense_variant", 
                          "protein_altering", "protein_coding"]

syn_coding_variants = ["synonymous", "stop_retained", "coding_sequence"]

noncoding_variants = ["intron", "downstream", "upstream", "non_coding", "intergenic"]

regulatory_variants = ["UTR", "regulatory_region", "mature_miRNA"]

splice_variants = ["splice"]

chrm_list = [str(i+1) for i in range(22)] + ["X", "Y"]

def vcf_read(filename, all_chrm = False, picklefile = None, skip = 0, pre_chrm = 0):
  try:
    with open(filename, 'r') as f:
      j = 0
      for _ in range(900+skip):
        if _ >= 900:
          j += 1
        next(f)
        if (j>0 and j%100000 == 0):
          print("%i skipped" %j)
      all_variants = [] if not all_chrm else {chrm: [] for chrm in chrm_list}
      i = skip
      prev_chrm = 0
      for line in f:
        tab_line = line.split('\t')
        chrm = tab_line[0]
        rsid = tab_line[2]
        ref = tab_line[3]
        alt = tab_line[4]
        filt_pass = tab_line[6]
        var_line = tab_line[7].split(';')
        var_dict = {var_i.split("=")[0]: var_i.split("=")[1] for var_i in var_line if "=" in var_i}
        
        ac = float(var_dict["AC"])
        ac_afr = float(var_dict["AC_afr"])
        ac_amr = float(var_dict["AC_amr"])
        ac_asj = float(var_dict["AC_asj"])
        ac_eas = float(var_dict["AC_eas"])
        ac_fin = float(var_dict["AC_fin"])
        ac_nfe = float(var_dict["AC_nfe"])
        ac_sas = float(var_dict["AC_sas"])

        an = float(var_dict["AN"])
        an_afr = float(var_dict["AN_afr"])
        an_amr = float(var_dict["AN_amr"])
        an_asj = float(var_dict["AN_asj"])
        an_eas = float(var_dict["AN_eas"])
        an_fin = float(var_dict["AN_fin"])
        an_nfe = float(var_dict["AN_nfe"])
        an_sas = float(var_dict["AN_sas"])
        
        '''
        hom = float(var_dict["nhomalt"])
        hom_afr = float(var_dict["nhomalt_afr"])
        hom_amr = float(var_dict["nhomalt_amr"])
        hom_asj = float(var_dict["nhomalt_asj"])
        hom_eas = float(var_dict["nhomalt_eas"])
        hom_fin = float(var_dict["nhomalt_fin"])
        hom_nfe = float(var_dict["nhomalt_nfe"])
        hom_sas = float(var_dict["nhomalt_sas"])
        
        af = float(var_line[2].split("=")[1])
        af_afr = float(var_dict["AF_afr"])
        af_amr = float(var_dict["AF_amr"])
        af_asj = float(var_dict["AF_asj"])
        af_eas = float(var_dict["AF_eas"])
        af_fin = float(var_dict["AF_fin"])
        af_nfe = float(var_dict["AF_nfe"])
        af_sas = float(var_dict["AF_sas"])
        '''

        variant_info = var_line[-1]
        variant = None
        for variant_type in variant_types:
          if variant_type in variant_info:
            variant = variant_type
        if variant == None:
          print(var_line)
          variant = input("Input variant type")

        this_variant = dict(chrm=chrm, rsid=rsid, ref=ref, alt=alt, filt_pass=filt_pass,
                            ac=ac, ac_afr=ac_afr, ac_amr=ac_amr, ac_asj=ac_asj, 
                            ac_eas=ac_eas, ac_fin=ac_fin, ac_nfe=ac_nfe, ac_sas=ac_sas,
                            an=an, an_afr=an_afr, an_amr=an_amr, an_asj=an_asj, 
                            an_eas=an_eas, an_fin=an_fin, an_nfe=an_nfe, an_sas=an_sas,
                            #af=af, af_afr=af_afr, af_amr=af_amr, af_asj=af_asj, 
                            #af_eas=af_eas, af_fin=af_fin, af_nfe=af_nfe, af_sas=af_sas, 
                            variant=variant)

        if all_chrm:
          all_variants[chrm].append(this_variant)
          if prev_chrm != chrm:
            if picklefile and (prev_chrm in chrm_list[pre_chrm:]):
              df = pd.DataFrame(all_variants[prev_chrm])
              pickle.dump(df, open(prev_chrm+picklefile,'wb'))
              print(prev_chrm + " done!")
            prev_chrm = chrm
        else:
          all_variants.append(this_variant)

        i += 1
        if (i%100000) == 0:
          print("%i variants" %i)

  except Exception as e:
    print(i)
    print()
    print(var_line)
    raise e

  if all_chrm:
    if picklefile:
      df = pd.DataFrame(all_variants[chrm_list[-1]])
      pickle.dump(df, open(chrm_list[-1]+picklefile,'wb'))
      print(prev_chrm + " done!")
  elif picklefile:
    df = pd.DataFrame(all_variants)
    pickle.dump(df, open(picklefile,'wb'))

  return all_variants


def vcf_rows(df):
  is_nonsyn_code = ((df.variant.isin(nonsyn_coding_variants) & (df.filt_pass == 'PASS')))
  nonsyn_code_rows = df[is_nonsyn_code]

  is_syn_code = ((df.variant.isin(syn_coding_variants) & (df.filt_pass == 'PASS')))
  syn_code_rows = df[is_syn_code]

  is_noncode = ((df.variant.isin(noncoding_variants) & (df.filt_pass == 'PASS')))
  noncode_rows = df[is_noncode]

  syn_rows = df[(is_syn_code | is_noncode)]

  is_reg = ((df.variant.isin(regulatory_variants) & (df.filt_pass == 'PASS')))
  reg_rows = df[is_reg]

  is_splice = ((df.variant.isin(splice_variants) & (df.filt_pass == 'PASS')))
  splice_rows = df[is_splice]

  print(len(df))
  print(len(syn_rows))
  print(len(nonsyn_code_rows))

  return syn_rows, nonsyn_code_rows

annotation = ['3_prime_UTR_variant',
              'stop_retained_variant',
              'missense_variant',
              'synonymous_variant',
              'frameshift_variant',
              'stop_gained',
              'inframe_deletion',
              'splice_region_variant',
              'intron_variant',
              'splice_donor_variant',
              'inframe_insertion',
              'splice_acceptor_variant',
              '5_prime_UTR_variant',
              'coding_sequence_variant',
              'stop_lost',
              'protein_altering_variant']

splice = ['splice_region_variant',
          'splice_donor_variant',
          'splice_acceptor_variant']

syn_code = ['synonymous_variant',
            'coding_sequence_variant',
            'stop_retained_variant']

utr = ['3_prime_UTR_variant',
       '5_prime_UTR_variant']

intron = ['intron_variant']

nonsyn_code = ['missense_variant',
               'frameshift_variant',
               'stop_gained',
               'inframe_deletion',
               'inframe_insertion',
               'protein_altering_variant',
               'stop_lost']


def csv_read_rows(filename):
  df = pd.read_csv(filename)

  is_splice = df.Annotation.isin(splice)
  splice_rows = df[is_splice]

  is_syn_code = df.Annotation.isin(syn_code)
  syn_code_rows = df[is_syn_code]

  is_utr = df.Annotation.isin(utr)
  utr_rows = df[is_utr]
  
  is_intron = df.Annotation.isin(intron)
  intron_rows = df[is_intron]
  
  is_nonsyn_code = df.Annotation.isin(nonsyn_code)
  nonsyn_code_rows = df[is_nonsyn_code]

  is_syn = ((is_syn_code | is_intron) | is_utr)
  syn_rows = df[is_syn]
  
  print(len(df))
  print(len(syn_rows))
  print(len(nonsyn_code_rows))

  return syn_rows, nonsyn_code_rows


def csv_graph(syn_rows, nonsyn_code_rows, genename, binsize):

  for pop in namepops:
    synx = np.log(syn_rows["Allele Count "+pop].divide(syn_rows["Allele Number "+pop]))
    nonsynx = np.log(nonsyn_code_rows["Allele Count "+pop].divide(nonsyn_code_rows["Allele Number "+pop]))
    #print(genename, pop, synx.count(), nonsynx.count())
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=binsize), opacity=.5))
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=binsize), opacity=.5))
    fig.update_layout(
            barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="P(x) (log scale))",
                type="log",
                range=[-14,0]
            ),
            xaxis=go.layout.XAxis(
                title_text="log(x)",
                range=[-14,0]
            ),
            title_text="%s bin=%1.2f (%s)" %(genename, binsize, pop)
        )
    fig.show()
    fig.write_image("%s bin=%1.2f (%s).png" %(genename, binsize, pop))

def vcf_graph(rows, pop, binsize, title, filename):
  expx = rows["ac_%s" %pop].divide(rows["an_%s" %pop])
  expx = expx[((expx > 0) & (expx < 1))]
  alphax, betax, _, _ = beta.fit(expx, floc=0, fscale=1)
  x = np.arange(0,1,binsize)
  binx = [(x[i+1] + x[i])/2 for i in range(len(x)-1)]
  y = [btdtr(alphax, betax, x[i+1]) - btdtr(alphax, betax, x[i]) for i in range(len(x)-1)]
  fig = go.Figure();
  fig.add_trace(go.Histogram(x=expx, histnorm='probability', name="Experimental", autobinx=False,
                             xbins=dict(start=0, end=1, size=binsize), opacity=.9))
  fig.add_trace(go.Bar(x=binx, y=y, name="Theory", opacity=.9))
  fig.update_layout(
          autosize=False,
          width=800,
          height=600,
          yaxis=go.layout.YAxis(
              title_text="P(x)",
              range=[0,1]
          ),
          xaxis=go.layout.XAxis(
              title_text="x",
              range=[0,1]
          ),
          title_text=title,
          legend_orientation="h"
      )
    #fig.write_image(filename)
  return (alphax, betax, kstest(expx, 'beta', args = (alphax, betax)).pvalue, expx.size)


def old_vcf_graph(syn_rows, nonsyn_code_rows, chrm, binsize):
  pop_vars = {}
  for i, pop in enumerate(pops):
    #synp = syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop])
    #synp1p = synp*(1-synp)
    #synx = np.log(synp1p)
    synx = syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop])
    synx = synx[synx > 0]
    synx = synx[synx < 1]
    #synx = np.log(syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop]))
    syna, synb, _, _ = beta.fit(synx, floc=0, fscale=1)
    x = np.arange(0,1,binsize)
    binx = [(x[i+1] + x[i])/2 for i in range(len(x)-1)]
    syny = [btdtr(syna, synb, x[i+1]) - btdtr(syna, synb, x[i]) for i in range(len(x)-1)]
    #syny =  beta.pdf(x, syna, synb)
    #print(chrm, namepops[i], synx.count(), nonsynx.count())
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=0, end=1, size=binsize), opacity=1))
    fig.add_trace(go.Bar(x=binx, y=syny, name="Synonymous Theory", opacity=1))
    #fig.add_trace(go.Scatter(x=x, y=syny, name="Synonymous Beta Fit", opacity=1))
    #fig.add_trace(go.Scatter(x=x, y=nonsyny, name="Nonsynonymous Beta Fit", opacity=1))
    
    fig.update_layout(
            #barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="P(x)",
                #type="log",
                range=[0,1]
            ),
            xaxis=go.layout.XAxis(
                title_text="x",
                #type="log",
                range=[0,1]
            ),
            title_text="Chromosome %s bin=%1.2f (%s)" %(chrm, binsize, namepops[i]),
            legend_orientation="h"
        )
    fig.write_image("images_chroms/%s/syn_chrm_%s_bin_%1.2f_%s.png" %(chrm, chrm, binsize, pop))
    fig = go.Figure();
    #nonsynx = np.log(nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop]))
    nonsynx = nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop])
    nonsynx = nonsynx[nonsynx > 0]
    nonsynx = nonsynx[nonsynx < 1]
    nonsyna, nonsynb, _, _ = beta.fit(nonsynx, floc=0, fscale=1)
    nonsyny = [btdtr(nonsyna, nonsynb, x[i+1]) - btdtr(nonsyna, nonsynb, x[i]) for i in range(len(x)-1)]
    #nonsyny =  beta.pdf(x, nonsyna, nonsynb)
    #nonsynp = nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop])
    #nonsynp1p = nonsynp*(1-nonsynp)
    #nonsynx = np.log(nonsynp1p)
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=0, end=1, size=binsize), opacity=1))
    fig.add_trace(go.Bar(x=binx, y=nonsyny, name="Nonsynonymous Theory", opacity=1))
    fig.update_layout(
            #barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="P(x)",
                #type="log",
                range=[0,1]
            ),
            xaxis=go.layout.XAxis(
                title_text="x",
                #type="log",
                range=[0,1]
            ),
            title_text="Chromosome %s bin=%1.2f (%s)" %(chrm, binsize, namepops[i]),
            legend_orientation="h"
        )
    fig.write_image("images_chroms/%s/nonsyn_chrm_%s_bin_%1.2f_%s.png" %(chrm, chrm, binsize, pop))
    pop_vars[namepops[i]] = (syna, synb, kstest(synx, 'beta', args = (syna, synb)), nonsyna, nonsynb, kstest(nonsynx, 'beta', args = (nonsyna, nonsynb)).pvalue)
  return pop_vars
    
if __name__ == "__main__":
  #df = vcf_read(r'E:\Fall2019\HST508\final_proj\gnomad.exomes.r2.1.1.sites.vcf\gnomad.exomes.r2.1.1.sites.vcf', 
  #              all_chrm = True, picklefile='_all.pickle', skip = 0, pre_chrm = 0)
  #df = vcf_read(r'E:\Fall2019\HST508\final_proj\gnomad.genomes.r2.1.1.sites.20.vcf', 
  #              all_chrm = False, picklefile='20_genome.pickle', skip = 0, pre_chrm = 0)
  for chrm in chrm_list:
    df = pickle.load(open(('fixed_pickles/%s_all.pickle' %chrm), 'rb'))
    save_vars = []
    print("Chromosome %s" %chrm)
    syn_rows, nonsyn_code_rows = vcf_rows(df)
    for binsize in [.01]:
      for i, pop in enumerate(pops):
        title = "Chromosome %s bin=%1.2f (%s)" %(chrm, binsize, namepops[i])
        for rows, category in [(syn_rows, "Synonymous"), (nonsyn_code_rows, "Nonsynonymous")]:
          filename = "images_chroms/%s/%s_chrm_%s_bin_%1.2f_%s.png" %(chrm, category, chrm, binsize, pop)
          save_vars.append([category, namepops[i], *vcf_graph(rows, pop, binsize, title, filename)])

    with open("fit_params_%s.txt" %chrm, "w") as f:
      f.write("Category, Population, Alpha, Beta, KS pvalue, Number\n")
      for line in save_vars:
        pop, cat, alp, bet, ksp, n = line
        f.write("%s, %s, %.5f, %.5f, %.9f, %i\n" %(pop, cat, alp, bet, ksp, n))

    #vcf_graph(syn_rows, nonsyn_code_rows, chrm, .01)
    '''
    with open("fit_params_%s.txt" %chrm, "w") as f:
      f.write("Category, Population, Alpha, Beta, KS pvalue\n")
      for pop in namepops:
        f.write("Synonymous, %s, %.5f, %.5f, %.9f\n" %(pop, save_vars[pop][0], save_vars[pop][1], save_vars[pop][2].pvalue))
        f.write("Nonsynonymous, %s, %.5f, %.5f, %.9f\n" %(pop, save_vars[pop][0], save_vars[pop][1], save_vars[pop][2].pvalue))
    '''

  '''
  lct_syn, lct_nonsyn = csv_read_rows("gnomAD_v2.1.1_ENSG00000115850_2019_12_13_19_29_29.csv")
  syne1_syn, syne1_nonsyn = csv_read_rows("gnomAD_v2.1.1_ENSG00000131018_2019_11_25_03_29_39.csv")

  for binsize in [.1,.25, .5, 1]:
    csv_graph(lct_syn, lct_nonsyn, "LCT", binsize)
    csv_graph(syne1_syn, syne1_nonsyn, "SYNE1", binsize)

  '''

