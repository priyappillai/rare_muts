import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pickle


pops = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas']
namepops = ["African", "Latino", "Ashkenazi Jewish", "East Asian", "European (Finnish)", "European (non-Finnish)", "South Asian"]
acpops = ['AC_%s=' % pop for pop in pops]
anpops = ['AN_%s=' % pop for pop in pops]
afpops = ['AF_%s=' % pop for pop in pops]

variant_types = ["splice", "start_lost", "stop_lost", "stop_gained", "frameshift", 
                 "inframe_insertion", "inframe_deletion", "missense_variant", "UTR", 
                 "synonymous", "stop_retained", "intron", "downstream", "upstream", 
                 "protein_altering", "non_coding", "intergenic", "regulatory_region", 
                 "incomplete_terminal_codon", "mature_miRNA"]

nonsyn_coding_variants = ["start_lost", "stop_lost", "stop_gained", "frameshift", 
                          "inframe_insertion", "inframe_deletion", "missense_variant", 
                          "protein_altering"]

syn_coding_variants = ["synonymous", "stop_retained"]

noncoding_variants = ["intron", "downstream", "upstream", "non_coding", "intergenic"]

regulatory_variants = ["UTR", "regulatory_region", "mature_miRNA"]

splice_variants = ["splice"]

def vcf_read(filename, picklefile = None):
  with open(filename, 'r') as f:
    for _ in range(900):
      f.readline()
    variants = []
    for line in f:
      tab_line = line.split('\t')
      chrm = tab_line[0]
      rsid = tab_line[2]
      ref = tab_line[3]
      alt = tab_line[4]
      filt_pass = tab_line[6]
      var_line = tab_line[7].split(';')

      ac = float(var_line[0].split("=")[1])
      ac_afr = float(var_line[78].split("=")[1])
      ac_amr = float(var_line[158].split("=")[1])
      ac_asj = float(var_line[586].split("=")[1])
      ac_eas = float(var_line[177].split("=")[1])
      ac_fin = float(var_line[526].split("=")[1])
      ac_nfe = float(var_line[494].split("=")[1])
      ac_sas = float(var_line[121].split("=")[1])

      an = float(var_line[1].split("=")[1])
      an_afr = float(var_line[79].split("=")[1])
      an_amr = float(var_line[159].split("=")[1])
      an_asj = float(var_line[587].split("=")[1])
      an_eas = float(var_line[178].split("=")[1])
      an_fin = float(var_line[527].split("=")[1])
      an_nfe = float(var_line[495].split("=")[1])
      an_sas = float(var_line[122].split("=")[1])

      af = float(var_line[2].split("=")[1])
      af_afr = float(var_line[80].split("=")[1])
      af_amr = float(var_line[160].split("=")[1])
      af_asj = float(var_line[588].split("=")[1])
      af_eas = float(var_line[179].split("=")[1])
      af_fin = float(var_line[528].split("=")[1])
      af_nfe = float(var_line[496].split("=")[1])
      af_sas = float(var_line[123].split("=")[1])

      variant_info = var_line[-1]
      variant = None
      for variant_type in variant_types:
        if variant_type in variant_info:
          variant = variant_type
      if variant == None:
        print(var_line)
        raise ValueError

      variants.append(dict(chrm=chrm, rsid=rsid, ref=ref, alt=alt, filt_pass=filt_pass,
                           ac=ac, ac_afr=ac_afr, ac_amr=ac_amr, ac_asj=ac_asj, 
                           ac_eas=af_eas, ac_fin=af_fin, ac_nfe=ac_nfe, ac_sas=ac_sas,
                           an=an, an_afr=an_afr, an_amr=an_amr, an_asj=an_asj, 
                           an_eas=an_eas, an_fin=an_fin, an_nfe=an_nfe, an_sas=an_sas,
                           af=af, af_afr=af_afr, af_amr=af_amr, af_asj=af_asj, 
                           af_eas=af_eas, af_fin=af_fin, af_nfe=af_nfe, af_sas=af_sas, 
                           variant=variant))

  df = pd.DataFrame(variants)
  if picklefile:
    pickle.dump(df, open(picklefile,'wb'))
  return df



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

  #print(len(df))
  #print(len(syn_rows))
  #print(len(nonsyn_code_rows))

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
               'coding_sequence_variant',
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
  
  #print(len(df))
  #print(len(syn_rows))
  #print(len(nonsyn_code_rows))

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


def vcf_graph(syn_rows, nonsyn_code_rows, chrm, binsize):
  for i, pop in enumerate(pops):
    synx = np.log(syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop]))
    nonsynx = np.log(nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop]))
    #print(chrm, namepops[i], synx.count(), nonsynx.count())
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=-10, end=0, size=binsize), opacity=.5))
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=-10, end=0, size=binsize), opacity=.5))
    fig.update_layout(
            barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="P(x) (log scale)",
                type="log",
                range=[-14,0]
            ),
            xaxis=go.layout.XAxis(
                title_text="log(x)",
                range=[-14,0]
            ),
            title_text="Chromosome %s bin=%1.2f (%s)" %(chrm, binsize, namepops[i])
        )
    fig.show()
    fig.write_image("Chromosome %s bin=%1.2f (%s).png" %(chrm, binsize, namepops[i]))

#df = vcf_read('gnomad.exomes.r2.1.1.sites.22.vcf', '22.pickle')
df = pickle.load(open('22.pickle', 'rb'))
syn_rows, nonsyn_code_rows = vcf_rows(df)
for binsize in [.1,.25, .5, 1]:
  vcf_graph(syn_rows, nonsyn_code_rows, '22', binsize)

lct_syn, lct_nonsyn = csv_read_rows("gnomAD_v2.1.1_ENSG00000115850_2019_12_13_19_29_29.csv")
syne1_syn, syne1_nonsyn = csv_read_rows("gnomAD_v2.1.1_ENSG00000131018_2019_11_25_03_29_39.csv")

for binsize in [.1,.25, .5, 1]:
  csv_graph(lct_syn, lct_nonsyn, "LCT", binsize)
  csv_graph(syne1_syn, syne1_nonsyn, "SYNE1", binsize)

