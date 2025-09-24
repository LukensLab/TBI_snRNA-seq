import os
import scanpy as sc
import scvi

input_file = '/path/to/data/soupx_combined.h5ad'
output_dir = '/path/to/data/solo_output'
os.makedirs(output_dir, exist_ok=True)

output_file_csv = 'solo_predictions.csv'
output_file_h5ad = os.path.join(output_dir, 'soupx_combined_with_solo.h5ad')
vae_dir = os.path.join(output_dir, 'vae_model')
solo_dir = os.path.join(output_dir, 'solo_model')
os.makedirs(vae_dir, exist_ok=True)
os.makedirs(solo_dir, exist_ok=True)

adata = sc.read_h5ad(input_file)
sc.pp.filter_genes(adata, min_cells=10)

if 'highly_variable' not in adata.var.columns:
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False, flavor='seurat_v3')

adata_hvg = adata[:, adata.var['highly_variable']].copy()

scvi.model.SCVI.setup_anndata(adata_hvg)
vae = scvi.model.SCVI(adata_hvg)
vae.train()
vae.save(vae_dir, overwrite=True)

solo = scvi.external.SOLO.from_scvi_model(vae, adata=adata_hvg, doublet_ratio=2)
solo.train()
solo.save(solo_dir, overwrite=True)

# Generate predictions and include sample_id
results = solo.predict()
results['prediction'] = solo.predict(soft=False)
results['sample_id'] = adata_hvg.obs['sample_id'].values

results.to_csv(os.path.join(output_dir, output_file_csv), index=True)

print("SOLO finished. Results saved to:")
print(f"- CSV: {os.path.join(output_dir, output_file_csv)}")
print(f"- AnnData: {output_file_h5ad}")
print(f"- VAE model: {vae_dir}")
print(f"- SOLO model: {solo_dir}")