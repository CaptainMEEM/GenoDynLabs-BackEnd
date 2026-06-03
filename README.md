# GenoDynLabs Backend — v2 (full overhaul)

Turns an uploaded 23andMe `raw.txt` into a ~100-page annotated genomic report
and emails it. This is a ground-up rewrite of the data pipeline plus a
non-blocking, encrypted, horizontally-scalable job architecture.

The report **layout is unchanged** (the part you liked): cover, summary band,
pastel topic banners, function-titled article headers, dense 5-column tables
with the "Your Genotype" cell highlighted. Everything that changed is in *what
feeds that layout*.

---

## 1. What was wrong, and what fixed it

**The notes were bad because the old precompute discarded ~90% of SNPedia.**
It kept only `gene/chr/pos/orientation` per SNP and `geno/mag/repute/summary`
per genotype. Measured against the real dumps, **63% of genotype summaries are
boilerplate** ("common in clinvar", "normal"), so ~2 of every 3 rows had a blank
or filler note. Meanwhile the SNP pages carry **6,400 GWAS associations** (trait
+ risk allele + odds ratio + study title) and `trait_df.csv` is the full GWAS
Catalog with population frequencies — none of it reached the note.

The fix is three new pieces feeding the old layout:

| File | Role |
|------|------|
| `services/precompute_reference.py` *(rewritten)* | Mines **all four** dumps into one rsid-keyed bundle: gene(s), orientation, GMAF, every SNP-page GWAS association, every GWAS-Catalog association (with population RAFs), per-genotype mag/repute/summary, and the LD map. |
| `services/notes.py` *(new)* | Layered note synthesis. For every variant it stacks: real genotype summary → SNP-page GWAS (dosage-aware, "1 copy of the G risk allele, ↑ OR 1.31") → GWAS-Catalog assoc. → curated gene-function baseline → cleaned page text. **Guarantees a non-blank note on every row.** |
| `services/gene_labels.py` *(new)* | ~120 curated gene→function labels and baseline sentences. This is the **FUT2→"Vitamin B12 Status"** fix: article titles describe what the gene *does*, not the raw symbol. |

**Measured result on the Mendel test genome:** note coverage went from ~37% to
**100%**, with ~5,500 notes citing a real odds ratio, and labels like
`Vitamin B12 Status (FUT2)`, `Folate Metabolism (MTHFR)`, `Alcohol Flush (ALDH2)`.

**Odds-ratio sanity:** values outside `[0.1, 10]` are effect sizes on continuous
traits (biomarker levels), not true case/control ORs, so the note shows a
direction arrow only — no misleading "OR 49.77".

### ~100 pages, stably

`report.build_pdf(..., target_pages=100)` uses an **article-first global budget**:
it ranks whole gene/trait sections by their strongest variant and fills up to a
row budget, collapsing the weak single-SNP tail into one dense "Additional"
table per topic. This lands at **94–99 pages whether the input is a sequenced
genome (12.8k matches) or a sparse consumer chip (5.7k matches)** — page count
no longer swings with input density. Tune with `target_pages` or
`REPORT_TARGET_PAGES`.

> **Honest caveat (from ModernPromethease's own author):** SNPedia per-SNP
> summaries are often outdated, and SNP-first reporting under-represents how
> polygenic most traits are. The report leans on GWAS-Catalog odds ratios where
> possible and carries a prominent disclaimer, but it remains SNPedia-grade data.

---

## 2. Architecture: non-blocking, concurrent, low Railway cost

```
                       (encrypted genome payload in Redis)
  browser ──TLS──> [ web: app.py ] ──enqueue──> Redis ──> [ worker: worker.py ] ──> email
                    auth, validate,             (RQ)       decrypt, annotate,
                    encrypt, return 202                    render PDF, send
```

* **Web never blocks.** `/upload` only authenticates, validates, encrypts, and
  enqueues — it returns `202 + job_id` in milliseconds. All heavy work runs in
  the worker. Poll `/status/<job_id>`; the email is the primary delivery.
* **Mass scale by replicas.** Workers share only the read-only reference bundle
  and have no mutable shared state, so concurrency = number of worker replicas.
  Scale web and worker independently in Railway.
* **Low resource use.**
  * Runtime image is **pandas/numpy-free** (the request path never parses a CSV).
  * The reference bundle (**6 MB**) loads **once per worker process** and is
    reused for every job; per-job cost is just the user's records.
  * Fully **in-memory** job: genome and PDF never touch disk — no temp files,
    no cleanup, no I/O contention.
  * Typical render ≈ **20–25 s** for ~100 pages in the background worker.

If `REDIS_URL` is unset, `app.py` falls back to an in-process thread pool for
local dev only (single web process; not for production scale).

---

## 3. Encryption — what it does and doesn't guarantee

A report can only be made by decrypting the genome and annotating it, so the
worker must see plaintext in memory; true zero-knowledge "server can never read
it" E2E is therefore impossible for the annotation step. What v2 *does*
guarantee (`services/crypto.py`):

* **In transit:** TLS everywhere (Railway HTTPS; Microsoft Graph HTTPS).
* **At rest in the queue:** the genome is wrapped in **authenticated AES-256-GCM
  before it ever enters Redis**. A Redis dump/backup/breach yields only
  ciphertext. Supports **key rotation** (list multiple keys; first encrypts, all
  are tried to decrypt).
* **In processing:** decrypt into memory only; **nothing is written to disk**, so
  there is no plaintext at rest and nothing to scrub.
* **Optional client-side sealed box:** the browser can fetch the server's X25519
  public key from `/pubkey` and **seal the genome to it before upload** (X25519 +
  HKDF + AES-256-GCM). Then the upload body is opaque to any proxy/CDN/request
  log and only this server's private key can open it. `app.py` auto-detects a
  sealed payload; plaintext `.txt` over TLS still works if you don't wire this up.

This is the strongest practical posture short of doing genomics in the browser.

---

## 4. Deploy on Railway

Two services from **one image** (Railway just overrides the start command):

| Service | Start command |
|---------|---------------|
| **web** | `gunicorn app:app --bind 0.0.0.0:$PORT --timeout 120 --workers 2 --threads 4` (default `CMD`) |
| **worker** | `python worker.py` |

Add the **Redis** plugin (sets `REDIS_URL`). Scale `worker` replicas for throughput.

### Environment variables

```
# Queue (Redis plugin sets this)
REDIS_URL=redis://...

# Encryption at rest (REQUIRED in prod). Generate with `make keys`.
GENODYN_DATA_KEYS=<base64 32-byte key>[,<older key for rotation>]

# Optional sealed-box upload path. Generate with `make keys`.
GENODYN_X25519_PRIV=<base64 X25519 private key>

# Auth + email (unchanged from v1)
FIREBASE_CREDENTIALS_JSON=<service account JSON or path>
AZURE_TENANT_ID=...
AZURE_CLIENT_ID=...
AZURE_CLIENT_SECRET=...
EMAIL_FROM=DONOTREPLY@genodynlabs.com

# Optional knobs
ALLOWED_ORIGINS=https://genodynlabs.com,https://genodynlabs.vercel.app
REPORT_TARGET_PAGES=100
WEB_WORKERS=2
WEB_THREADS=4
```

Generate keys locally:

```
make keys
# prints GENODYN_DATA_KEYS=... and the X25519 private/public pair
```

### API

* `POST /upload` (auth) — multipart `file`; plaintext `.txt` or a sealed blob.
  → `202 {job_id, status, emailed_to}`
* `GET  /status/<job_id>` (auth) → `{status: queued|processing|done|failed}`
* `GET  /pubkey` → `{x25519_public_key}` for the optional client seal
* `GET  /health` → `{status, dispatch, encryption}`

---

## 5. Regenerating the reference data

The deploy ships the prebuilt `data/reference.pkl` (6 MB); the source CSVs are
git-ignored (they're only needed to rebuild). When you refresh the SNPedia/GWAS
dumps, drop the four CSVs into `data/` and:

```
make reference        # parses CSVs -> data/reference.pkl  (~10 s)
git add data/reference.pkl && git commit -m "refresh reference data"
```

CSV layout expected in `data/`: `snp_df.csv`, `geno_df.csv`, `trait_df.csv`,
`equilibrium_df.csv` (the ModernPromethease formats).

To grow coverage of curated labels/notes, just add entries to
`GENE_FUNCTION` in `services/gene_labels.py` — no data rebuild needed.

---

## 6. File map

```
app.py                          web service (cheap, non-blocking path only)
worker.py                       RQ worker (warms bundle, serves the queue)
services/
  precompute_reference.py       REWRITTEN — mines all 4 dumps -> reference.pkl
  annotate.py                   REWRITTEN — in-memory, rich records
  notes.py                      NEW       — layered note synthesis (100% coverage)
  gene_labels.py                NEW       — curated function labels + baselines
  report.py                     same layout, rich data + page budget
  topics.py                     topic taxonomy (extended)
  crypto.py                     NEW       — AES-256-GCM at rest + X25519 sealed box
  jobs.py                       REWRITTEN — decrypt->annotate->render->email, in-mem
  firebase_auth.py              unchanged
  email_sender.py               unchanged
data/reference.pkl              prebuilt bundle (ships in image)
requirements.txt                runtime deps (no pandas)
requirements-build.txt          build-only deps (pandas) for `make reference`
Dockerfile / Procfile / Makefile
```
