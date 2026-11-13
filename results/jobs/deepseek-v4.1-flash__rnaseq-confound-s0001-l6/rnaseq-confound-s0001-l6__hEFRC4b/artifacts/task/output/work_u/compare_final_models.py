import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts=pd.read_csv('/task/output/work_u/gene_counts_merged.tsv',sep='\t',index_col=0)
sheet=pd.read_csv('/task/input/samples.tsv',sep='\t',dtype=str).set_index('sample_id')
cond=sheet.condition.copy(); cond['RS26']='treated'; cond['RS75']='vehicle'
sex={'RS22':'M','RS26':'F','RS27':'F','RS47':'M','RS62':'M','RS69':'F','RS75':'M','RS79':'M','RS80':'F','RS89':'F','RS97':'F'}
actual_batch=sheet.library_prep_batch.copy(); actual_batch['RS26']=sheet.loc['RS75','library_prep_batch']; actual_batch['RS75']=sheet.loc['RS26','library_prep_batch']
use=[s for s in sheet.index if s!='RS19']
mat=counts[use].T.round().astype(int)
mail=pd.read_csv('/task/input/previous_analysis/mail_top_responders.tsv',sep='\t',dtype=str)

for name,design in [('sexonly','~condition + sex'),('actualbatchsex','~condition + session + sex'),('sheetbatchsex','~condition + session + sex')]:
    meta=pd.DataFrame({'condition':cond.loc[use].values,'sex':[sex[s] for s in use]},index=use)
    if name=='actualbatchsex':
        meta['session']=['s'+actual_batch[s].replace('-','') for s in use]
    elif name=='sheetbatchsex':
        meta['session']=['s'+sheet.loc[s,'library_prep_batch'].replace('-','') for s in use]
    dds=DeseqDataSet(counts=mat,metadata=meta,design=design,quiet=True,n_cpus=2)
    dds.deseq2()
    st=DeseqStats(dds,contrast=['condition','treated','vehicle'],alpha=0.05,quiet=True,n_cpus=2)
    st.summary()
    res=st.results_df.reindex(counts.index)
    out=pd.DataFrame({'gene_id':counts.index,'log2fc':res['log2FoldChange'].values,'pvalue':res['pvalue'].values,'padj':res['padj'].values})
    out.to_csv(f'/task/output/work_u/de_{name}.tsv',sep='\t',index=False,na_rep='NA')
    sig=set(out.loc[out.padj<0.05,'gene_id'])
    print(name,'tested',out.padj.notna().sum(),'sig',len(sig),'mail confirmed',mail.gene_id.isin(sig).sum(),'/',len(mail))
    print(out.sort_values('padj').head(15).to_string(index=False))
