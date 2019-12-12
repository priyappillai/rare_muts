import pandas as pd
import numpy as np
import plotly.graph_objects as go


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

with open('gnomad.exomes.r2.1.1.sites.22.vcf', 'r') as f:
  for _ in range(900):
    f.readline()
  variants = []
  for line in f:
    tab_line = line.split('\t')
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

    variants.append(dict(rsid=rsid, ref=ref, alt=alt, filt_pass=filt_pass,
                         ac=ac, ac_afr=ac_afr, ac_amr=ac_amr, ac_asj=ac_asj, 
                         ac_eas=af_eas, ac_fin=af_fin, ac_nfe=ac_nfe, ac_sas=ac_sas,
                         an=an, an_afr=an_afr, an_amr=an_amr, an_asj=an_asj, 
                         an_eas=an_eas, an_fin=an_fin, an_nfe=an_nfe, an_sas=an_sas,
                         af=af, af_afr=af_afr, af_amr=af_amr, af_asj=af_asj, 
                         af_eas=af_eas, af_fin=af_fin, af_nfe=af_nfe, af_sas=af_sas, 
                         variant=variant))

df = pd.DataFrame(variants)
print(df)

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

for i, pop in enumerate(pops):
    synx = syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop]) 
    #np.log(syn_rows["ac_%s" %pop].divide(syn_rows["an_%s" %pop]))
    nonsynx = nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop]) 
    #np.log(nonsyn_code_rows["ac_%s" %pop].divide(nonsyn_code_rows["an_%s" %pop]))
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=0, end=1, size=.005), opacity=.75))
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=0, end=1, size=.005), opacity=.75))
    fig.update_layout(
            barmode='overlay',
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
                range=[0,1]
            ),
            title_text="%s - Synonymous vs Nonsynonymous" %namepops[i]
        )
    fig.show()
    fig.write_image("%s_non_syn.png" %pop)

