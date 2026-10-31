#!/usr/bin/env python3
"""Reaggregate verified GSE311899 epithelial-origin 10x matrices by patient and stage.

The common disjoint epithelial identity gate is deliberately applied to raw
counts. This is a paired longitudinal extension, not an independent replicate.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, re, sys, tarfile
from pathlib import Path
import numpy as np
from scipy import io, sparse

EPITHELIAL=("EPCAM","KRT8","KRT18","KRT19","KRT7","MUC1","PAX8","FOXJ1","PIFO","TPPP3","CAPS")
NON_EPITHELIAL=("DCN","LUM","COL1A1","COL3A1","PTPRC","LST1","PECAM1","VWF","RGS5","TAGLN","MYH11")
MIN_MARKERS, MIN_IDENTITY, MIN_CELLS=3,0.5,20
def sha(p:Path):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
 return h.hexdigest()
def member(tar,suffix):
 x=[m for m in tar.getmembers() if m.name.endswith('/'+suffix)]
 if len(x)!=1:raise ValueError(f'{suffix}: expected 1 member, found {len(x)}')
 return x[0]
def readlines(tar,m):
 f=tar.extractfile(m)
 if f is None:raise ValueError('unreadable tar member')
 with gzip.GzipFile(fileobj=f,mode='rb') as g:return [z.decode('utf-8',errors='replace').strip() for z in g]
def load_10x(path:Path):
 with tarfile.open(path,'r:gz') as tar:
  fm=member(tar,'features.tsv.gz'); bm=member(tar,'barcodes.tsv.gz'); mm=member(tar,'matrix.mtx.gz')
  features=readlines(tar,fm); barcodes=readlines(tar,bm); f=tar.extractfile(mm)
  if f is None:raise ValueError('matrix member unreadable')
  with gzip.GzipFile(fileobj=f,mode='rb') as g: mat=io.mmread(g).tocsr()
 if mat.shape != (len(features),len(barcodes)):raise ValueError(f'10x dimensions do not match member lines: {mat.shape}, {len(features)}, {len(barcodes)}')
 if mat.nnz and (mat.data.min()<0 or not np.allclose(mat.data,np.rint(mat.data),atol=1e-8,rtol=0)):raise ValueError('matrix has non-integer/negative values')
 symbols=np.array([x.split('\t')[1] if len(x.split('\t'))>=2 else '' for x in features],dtype=object)
 if np.any(symbols==''):raise ValueError('feature table lacks symbol values')
 unique, inv=np.unique(symbols,return_inverse=True)
 if len(unique)!=len(symbols):
  grouping=sparse.csr_matrix((np.ones(len(symbols)),(inv,np.arange(len(symbols)))),shape=(len(unique),len(symbols)))
  mat=(grouping@mat).tocsr()
 return unique,mat,len(barcodes),len(features)
def stage_short(s):return 'pre' if s=='Pre-AS therapy' else ('post' if s=='Post-AS therapy' else 'invalid')
def main(root:Path,source:Path):
 recpath=root/'metadata'/'external'/'GSE311899_epithelial_archive_receipt.tsv'; outcounts=root/'interim'/'GSE311899_paired_disjoint_epithelial_pseudobulk_counts.tsv.gz'; outmeta=root/'metadata'/'GSE311899_paired_disjoint_epithelial_pseudobulk_metadata.tsv'; outoverlap=root/'results'/'GSE311899_GSE215968_identifier_overlap_audit.tsv'; outqa=root/'metadata'/'GSE311899_paired_disjoint_epithelial_qa.json'
 if any(p.exists() for p in(outcounts,outmeta,outoverlap,outqa)):raise FileExistsError('Refusing to overwrite P2 outputs')
 rows=list(csv.DictReader(recpath.open(encoding='utf-8'),delimiter='\t'))
 ledger=list(csv.DictReader((root/'metadata'/'external'/'GSE311899_sample_ledger.tsv').open(encoding='utf-8'),delimiter='\t'))
 library_titles={r['geo_sample']:r['library_title'] for r in ledger}
 paired_ids=sorted({r['individual'] for r in rows if r['treatment_stage']=='Pre-AS therapy'} & {r['individual'] for r in rows if r['treatment_stage']=='Post-AS therapy'})
 rows=[r for r in rows if r['individual'] in paired_ids]
 if len(rows)!=2*len(paired_ids):raise ValueError('Paired archive ledger is not one pre + one post per patient')
 all_pb=[]; met=[]; failed_libraries=[]; symbols0=None
 for r in sorted(rows,key=lambda x:(x['individual'],x['treatment_stage'])):
  path=root/r['expected_local_path']; symbols,mat,barcode_n,feature_n=load_10x(path)
  if symbols0 is None:symbols0=symbols
  elif not np.array_equal(symbols0,symbols):raise ValueError(f'Feature symbols differ across archives: {path.name}')
  ix={g:i for i,g in enumerate(symbols)}; need=set(EPITHELIAL+NON_EPITHELIAL)
  if need-set(ix):raise ValueError(f'Missing identity genes in {path.name}: {sorted(need-set(ix))}')
  libs=np.asarray(mat.sum(axis=0)).ravel(); ep=mat[[ix[g] for g in EPITHELIAL],:].toarray(); ne=mat[[ix[g] for g in NON_EPITHELIAL],:].toarray()
  if np.any(libs<=0):raise ValueError(f'Zero-library filtered cell in {path.name}')
  detected=(ep>=1).sum(axis=0); ident=np.log1p(ep/libs*10000).mean(axis=0)-np.log1p(ne/libs*10000).mean(axis=0); passed=(detected>=MIN_MARKERS)&(ident>=MIN_IDENTITY)
  patient=re.search(r'(\d+)$',r['individual']).group(1); stage=stage_short(r['treatment_stage'])
  record={'pseudobulk_id':f'P{patient}_{stage}','patient_id':f'Patient {patient}','patient_token':patient,'stage':stage,'geo_sample':r['geo_sample'],'library_title':library_titles[r['geo_sample']],'source_url':r['source_url'],'raw_archive_sha256':r['sha256'],'filtered_matrix_barcodes':barcode_n,'features_before_symbol_deduplication':feature_n,'features_after_symbol_deduplication':len(symbols),'identity_gate_candidate_cells':barcode_n,'identity_gate_retained_cells':int(passed.sum()),'identity_gate_pass_fraction':float(passed.mean()),'median_identity_difference_mean_logCP10k':float(np.median(ident)),'matrix':'10x filtered_feature_bc_matrix raw integer sum'}
  if int(passed.sum()) < MIN_CELLS:
   failed_libraries.append(record)
   continue
  all_pb.append(np.asarray(mat[:,passed].sum(axis=1)).ravel().astype(np.int64)); met.append(record)
 retained_patients={x['patient_id'] for x in met if x['stage']=='pre'} & {x['patient_id'] for x in met if x['stage']=='post'}
 kept=[i for i,x in enumerate(met) if x['patient_id'] in retained_patients]
 met=[met[i] for i in kept]; all_pb=[all_pb[i] for i in kept]
 if len(retained_patients)<3: raise ValueError(f'Fewer than three paired patients retained after fixed {MIN_CELLS}-cell gate: {len(retained_patients)}')
 order=np.argsort([x['pseudobulk_id'] for x in met]); met=[met[i] for i in order]; pbs=[all_pb[i] for i in order]; ids=[x['pseudobulk_id'] for x in met];
 with gzip.open(outcounts,'wt',encoding='utf-8') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['gene']+ids)
  for i,g in enumerate(symbols0):w.writerow([g]+[int(x[i]) for x in pbs])
 with outmeta.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(met[0]),delimiter='\t');w.writeheader();w.writerows(met)
 # De-duplicate role by identifier token; source IDs are retained rather than assuming exact library identity.
 origin=[]
 for r in csv.DictReader((source/'results'/'GSE215968_origin_sample_counts.tsv').open(encoding='utf-8'),delimiter='\t'):
  if r.get('cohort_role')=='primary_local' and r.get('condition')=='AS':origin.append(r['donor_id'])
 eligible={r['donor_id'] for r in csv.DictReader((root/'metadata'/'GSE215968_disjoint_identity_gate_pseudobulk_metadata.tsv').open(encoding='utf-8'),delimiter='\t') if r.get('condition')=='AS'}
 outrows=[]
 for p in paired_ids:
  token=int(re.search(r"(\d+)$",p).group(1))
  primary=f'AS_{token:02d}'
  outrows.append({'GSE311899_patient':p,'eligible_for_fixed_gate_paired_analysis':str(p in retained_patients),'GSE215968_identifier_matched_AS_donor':primary,'in_GSE215968_primary_AS_ledger':str(primary in origin),'in_GSE215968_disjoint_primary_analysis':str(primary in eligible),'overlap_decision':'not_an_independent_replication; identifier-token overlap treated conservatively as overlapping patients','exact_library_identity':'not_established_by_identifier_only'})
 with outoverlap.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter='\t');w.writeheader();w.writerows(outrows)
 qa={'analysis_id':'IUA_RJ_P2_LONGITUDINAL','paired_patient_n_before_fixed_gate':len(paired_ids),'paired_patient_n_after_fixed_gate':len(retained_patients),'patients_after_fixed_gate':sorted(retained_patients),'failed_libraries_below_fixed_gate':failed_libraries,'raw_counts_nonnegative_integer':True,'gene_symbol_deduplication':'sum raw counts of repeated symbols within each 10x matrix before pseudobulk','identity_gate':{'epithelial_markers':list(EPITHELIAL),'disjoint_non_epithelial_markers':list(NON_EPITHELIAL),'minimum_detected_markers':MIN_MARKERS,'minimum_identity_difference_mean_logCP10k':MIN_IDENTITY,'minimum_retained_cells_per_library':MIN_CELLS},'overlap_status':'Patient identifier tokens overlap GSE215968 AS donors; P2 is a longitudinal extension and must not be counted as independent disease replication. Exact library identity was not inferred from identifiers alone.','nonclaims':['no treatment efficacy estimate','no causal effect','no clinical endpoint linkage','no independent replication']}
 outqa.write_text(json.dumps(qa,indent=2)+"\n",encoding='utf-8');print(json.dumps(qa,indent=2))
if __name__=='__main__':
 if len(sys.argv)!=3:raise SystemExit('Usage: 12_reaggregate_gse311899_paired_epithelial.py REVISION_ROOT SOURCE_ROOT')
 main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
