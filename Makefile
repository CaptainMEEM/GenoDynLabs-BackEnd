# Regenerate the reference bundle from the SNPedia/GWAS CSVs in data/ and,
# optionally, the deana evidence-pack shards for nutrient-level enrichment.
# Run this whenever you refresh source data, then commit data/reference.pkl.
#
# Set DEANA_SHARDS to the path of deana's evidence-pack shards dir to fold in
# nutrient level/status associations (only nutrient-relevant records are kept,
# so the bundle stays small). Leave unset to build from SNPedia data alone.
DEANA_SHARDS ?=
DEANA_LEVELS ?= high,moderate

.PHONY: reference keys
reference:
	pip install -r requirements-build.txt
	python -m services.precompute_reference \
	  data/snp_df.csv data/geno_df.csv data/trait_df.csv \
	  data/equilibrium_df.csv data/reference.pkl $(DEANA_SHARDS) $(DEANA_LEVELS)

# Print a fresh AES-256 data key and an X25519 keypair for the env vars.
keys:
	@echo "GENODYN_DATA_KEYS=$$(python -m services.crypto keygen)"
	@python -m services.crypto x25519
