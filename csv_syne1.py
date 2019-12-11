import pandas as pd
import numpy as np
import plotly.graph_objects as go
import allel


callset = allel.read_vcf('gnomad.exomes.r2.1.1.sites.22.vcf', fields=['variants/vep'])#, 'variants/variant_type', 'variants/allele_type', 'variants/AC', 'variants/AF'])

df = pd.read_csv("gnomAD_v2.1.1_ENSG00000131018_2019_11_25_03_29_39.csv")

pops = ["African", "Latino", "Ashkenazi Jewish", "East Asian", "European (Finnish)", "European (non-Finnish)", "Other", "South Asian"]

annotations = ['3_prime_UTR_variant', 'stop_retained_variant',
       'synonymous_variant', 'missense_variant', 'inframe_insertion',
       'stop_gained', 'inframe_deletion', 'frameshift_variant',
       'splice_acceptor_variant', 'intron_variant',
       'splice_region_variant', 'splice_donor_variant', 'start_lost',
       '5_prime_UTR_variant', 'coding_sequence_variant', 'stop_lost',
       'protein_altering_variant']

#Alt: check if there is a protein consequence/if consequence starts with "c" 
coding = ['stop_retained_variant', 'synonymous_variant',
          'missense_variant', 'inframe_insertion', 'stop_gained', 
          'inframe_deletion', 'frameshift_variant', 'start_lost', 
          'coding_sequence_variant', 'stop_lost', 
          'protein_altering_variant']

noncoding = ['3_prime_UTR_variant','splice_acceptor_variant', 'intron_variant',
             'splice_donor_variant','5_prime_UTR_variant']

#to check if coding, check if there is a protein consequence (will have "c." and "(p.=)")
#probably synonymous but not clear, there are some LoFs??
confusing = 'splice_region_variant'

#Alt: if coding, check if protein consequence has "(p.=)"/starts with "c"
synonymous = ['3_prime_UTR_variant', 'stop_retained_variant',
              'synonymous_variant', 'intron_variant', '5_prime_UTR_variant']

nonsynonymous = ['missense_variant', 'inframe_insertion', 'stop_gained',
                 'inframe_deletion', 'frameshift_variant',
                 'splice_acceptor_variant', 'splice_donor_variant',
                 'start_lost', 'coding_sequence_variant', 'stop_lost',
                 'protein_altering_variant']


syn_code = ['stop_retained_variant', 'synonymous_variant']
is_syn_code = (df.Annotation.isin(syn_code) | \
              ((df.Annotation == 'splice_region_variant') & \
               (df["Protein Consequence"].notna())))
syn_code_rows = df[is_syn_code]

nonsyn_code = ['missense_variant', 'inframe_insertion', 'stop_gained',
               'inframe_deletion', 'frameshift_variant', 'start_lost', 
               'coding_sequence_variant', 'stop_lost', 'protein_altering_variant']
is_nonsyn_code = (df.Annotation.isin(nonsyn_code))
nonsyn_code_rows = df[is_nonsyn_code]

syn_noncode = ['3_prime_UTR_variant', 'intron_variant', '5_prime_UTR_variant']
is_syn_noncode = (df.Annotation.isin(syn_noncode) | \
                 ((df.Annotation == 'splice_region_variant') & \
                  (df["Protein Consequence"].isna())))
syn_noncode_rows = df[is_syn_noncode]

nonsyn_noncode = ['splice_acceptor_variant', 'splice_donor_variant']
is_nonsyn_noncode = (df.Annotation.isin(nonsyn_noncode))
nonsyn_noncode_rows = df[is_nonsyn_noncode]

for pop in pops:
    synx = np.log(syn_code_rows["Allele Count "+pop].divide(syn_code_rows["Allele Number "+pop]))
    nonsynx = np.log(nonsyn_code_rows["Allele Count "+pop].divide(nonsyn_code_rows["Allele Number "+pop]))
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=.25), opacity=.75))
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=.25), opacity=.75))
    fig.update_layout(
            barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="log(P(x))",
                type="log",
                range=[-3,0]
            ),
            xaxis=go.layout.XAxis(
                title_text="log(x)",
                range=[-14,0]
            ),
            title_text="Coding ("+pop+")"
        )
    fig.show()
    fig.write_image("Coding ("+pop+")"+".png")

for pop in pops:
    synx = np.log(syn_noncode_rows["Allele Count "+pop].divide(syn_noncode_rows["Allele Number "+pop]))
    nonsynx = np.log(nonsyn_noncode_rows["Allele Count "+pop].divide(nonsyn_noncode_rows["Allele Number "+pop]))
    fig = go.Figure();
    fig.add_trace(go.Histogram(x=synx, histnorm='probability', name="Synonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=.25), opacity=.75))
    fig.add_trace(go.Histogram(x=nonsynx, histnorm='probability', name="Nonsynonymous", autobinx=False,
                               xbins=dict(start=-14, end=0, size=.25), opacity=.75))
    fig.update_layout(
            barmode='overlay',
            autosize=False,
            width=800,
            height=600,
            yaxis=go.layout.YAxis(
                title_text="P(x) (log scale)",
                type="log",
                range=[-3,0]
            ),
            xaxis=go.layout.XAxis(
                title_text="log(x)",
                range=[-14,0]
            ),
            title_text="Noncoding ("+pop+")"
        )
    fig.show()
    fig.write_image("Noncoding ("+pop+")"+".png")
