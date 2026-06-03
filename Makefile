# Regenerate the reference bundle from the SNPedia/GWAS CSVs in data/.
# Run this whenever you refresh the source data, then commit data/reference.pkl.
.PHONY: reference keys
reference:
	pip install -r requirements-build.txt
	python -m services.precompute_reference \
	  data/snp_df.csv data/geno_df.csv data/trait_df.csv \
	  data/equilibrium_df.csv data/reference.pkl

# Print a fresh AES-256 data key and an X25519 keypair for the env vars.
keys:
	@echo "GENODYN_DATA_KEYS=$$(python -m services.crypto keygen)"
	@python -m services.crypto x25519
