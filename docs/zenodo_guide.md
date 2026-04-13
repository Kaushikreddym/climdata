# Zenodo Archival & DOI Guide

This document provides instructions for archiving **climdata** on **Zenodo** to obtain a persistent Digital Object Identifier (DOI) for research reproducibility and citation purposes.

## Why Zenodo?

**Zenodo** is an open-access repository funded by CERN and the European Commission. Archiving climdata on Zenodo:

✅ **Creates a DOI** — Persistent identifier for stable citation
✅ **Enables reproducibility** — Exact versions are archived indefinitely
✅ **Improves discoverability** — Listed in CERN, OpenDOAR, and global registries
✅ **Supports open science** — Free access, no publication fees
✅ **Integrates with GitHub** — Automatic archival on each release

---

## Step 1: Create a Zenodo Account

1. Go to [https://zenodo.org/](https://zenodo.org/)
2. Click **"Sign up"** or login with GitHub
3. Complete your profile (name, email, affiliation)
4. Verify your email address

---

## Step 2: Enable GitHub Integration

### Option A: Automatic Archival on Release (Recommended)

1. Go to [Zenodo GitHub Settings](https://zenodo.org/account/settings/github/)
2. Click **"Connect"** to link your GitHub account
3. Grant Zenodo permission to access your repositories
4. Select the **climdata** repository from the list
5. Enable the toggle to activate automatic archival on releases

Now, every time you create a GitHub release, Zenodo will automatically:
- Download the release assets
- Create a new DOI
- Archive the version indefinitely

### Option B: Manual Upload

If you prefer manual control:

1. Go to [Zenodo Upload](https://zenodo.org/deposit/)
2. Click **"New Upload"**
3. Follow the instructions below (Step 3)

---

## Step 3: Prepare Release Metadata

Before creating a GitHub release, ensure your project has proper metadata files:

### ✅ CITATION.cff (Already Created)
- Located: `/CITATION.cff`
- Zenodo automatically reads this for citation information
- GitHub displays a "Cite this repository" button using this file

### ✅ codemeta.json (Already Created)
- Located: `/codemeta.json`
- Machine-readable metadata following CodeMeta standard
- Used by registries and search engines

### Update These Fields in CITATION.cff:

```yaml
version: 0.5.0  # Match your release tag
date-released: "2024-01-01"  # Update with release date
authors:
  - family-names: Muduchuru
    given-names: Kaushik
    orcid: "https://orcid.org/XXXX-XXXX-XXXX-XXXX"  # Add your ORCID
```

---

## Step 4: Create a GitHub Release

1. Go to [climdata Releases](https://github.com/Kaushikreddym/climdata/releases)
2. Click **"Create a new release"**
3. Enter:
   - **Tag version**: `v0.5.0` (must match pyproject.toml version)
   - **Release title**: `climdata v0.5.0`
   - **Description**: Include key features, bug fixes, and links

### Release Description Template:

```markdown
# climdata v0.5.0

## ✨ New Features
- Multi-provider weather data extraction (MSWX, CMIP, ERA5, DWD, HYRAS, POWER)
- Climate extreme indices calculation
- Bias correction and statistical downscaling (BCSD)
- Imputation for missing data gaps

## 🐛 Bug Fixes
- Fixed Google Drive API authentication
- Improved error handling for failed downloads

## 📚 Documentation
- Added MSWX dataset guide with licensing information
- Enhanced installation instructions
- New example notebooks

## 🔗 Links
- Documentation: https://Kaushikreddym.github.io/climdata
- MSWX Access: https://www.gloh2o.org/mswx/
- Changelog: https://github.com/Kaushikreddym/climdata/releases

### Citation

If you use climdata in your research, please cite:

```bibtex
@software{muduchuru2024climdata,
  title={climdata: Automated Climate Data Extraction and Processing},
  author={Muduchuru, Kaushik},
  year={2024},
  version={0.5.0},
  url={https://github.com/Kaushikreddym/climdata},
  doi={<DOI will be provided by Zenodo>}
}
```
```

4. Attach release artifacts (optional):
   - Source code ZIP/tarball
   - Wheels or distributions

5. Click **"Publish release"**

---

## Step 5: Verify Zenodo Archival

After creating a GitHub release (if auto-integration is enabled):

1. Go to [https://zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/)
2. Look for **climdata** in the list
3. Click on it to see archived versions
4. Each version will have a unique DOI

### Manual Verification:

1. Search [Zenodo](https://zenodo.org/search?q=climdata) for "climdata"
2. Your project should appear with all versions listed
3. Each version shows its DOI

---

## Step 6: Update Project with DOI

Once you have a DOI from Zenodo, update your project:

### Update README.md:

```markdown
[![DOI](https://zenodo.org/badge/XXXXX/climdata.svg)](https://zenodo.org/record/XXXXX)
```

Replace `XXXXX` with your Zenodo record ID (found in the DOI: `https://doi.org/10.5281/zenodo.XXXXX`).

### Update pyproject.toml:

```toml
[project.urls]
"DOI" = "https://doi.org/10.5281/zenodo.XXXXX"
"Zenodo" = "https://zenodo.org/record/XXXXX"
```

### Update CITATION.cff:

```yaml
identifiers:
  - description: "DOI for zenodo archive"
    type: doi
    value: "10.5281/zenodo.XXXXX"
```

---

## Step 7: Update Documentation Links

### In README.md:

Add a "Citation" section:

```markdown
## Citation

If you use climdata in your research, please cite:

**BibTeX:**
```bibtex
@software{muduchuru2024climdata,
  title={climdata: Automated Climate Data Extraction and Processing},
  author={Muduchuru, Kaushik},
  year={2024},
  version={0.5.0},
  url={https://github.com/Kaushikreddym/climdata},
  doi={10.5281/zenodo.XXXXX}
}
```

**APA:**
Muduchuru, K. (2024). climdata: Automated Climate Data Extraction and Processing (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXX

**BibLaTeX:**
```bibtex
@software{muduchuru2024climdata,
  title={climdata},
  author={Muduchuru, Kaushik},
  date={2024},
  version={0.5.0},
  url={https://github.com/Kaushikreddym/climdata},
  doi={10.5281/zenodo.XXXXX},
}
```
```

---

## Zenodo DOI Format

Your DOI will look like:
```
https://doi.org/10.5281/zenodo.XXXXXXX
```

- `10.5281` = Zenodo's DOI prefix (assigned by CERN)
- `zenodo.XXXXXXX` = Your unique record identifier

---

## Benefits of DOI

✅ **Permanent link** — Never broken, even if GitHub URL changes
✅ **Citable** — Use in papers, proposals, and CV
✅ **Discoverable** — Listed in Google Scholar, Crossref, DataCite
✅ **Version-specific** — Each release has its own DOI
✅ **Machine-readable** — Metadata harvestable by registries

---

## Publishing on PyPI with Zenodo

When you publish to PyPI, add Zenodo metadata:

```toml
[project]
name = "climdata"
version = "0.5.0"
# ... other fields ...

[project.urls]
"Homepage" = "https://github.com/Kaushikreddym/climdata"
"DOI" = "https://doi.org/10.5281/zenodo.XXXXX"
"Zenodo" = "https://zenodo.org/record/XXXXX"
```

---

## Verifying Your Zenodo Record

Check these registries to verify your software is discoverable:

1. **OpenDOAR** (https://www.opendoar.org/) — Directory of open-access repositories
2. **Zenodo** (https://zenodo.org/search?q=climdata)
3. **DataCite** (https://search.datacite.org/) — Global data citation index
4. **Software Heritage** (https://www.softwareheritage.org/) — Code preservation
5. **Google Scholar** (https://scholar.google.com/) — Academic search

---

## Using the DOI in Publications

### Journal Submission:

Include in "Data and Code Availability" section:

> The climdata software is available on GitHub (https://github.com/Kaushikreddym/climdata) and archived on Zenodo (https://doi.org/10.5281/zenodo.XXXXX).

### Grant Proposals:

> We will use climdata (Muduchuru, 2024; https://doi.org/10.5281/zenodo.XXXXX) for automated weather data extraction...

### PhD Thesis:

Include in bibliography with DOI:

```
Muduchuru, K. (2024). climdata: Automated Climate Data Extraction and Processing (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXX
```

---

## Troubleshooting

### Problem: GitHub release not appearing on Zenodo

**Solution:**
1. Ensure Zenodo GitHub integration is enabled in account settings
2. Use a semantic version tag (v0.5.0, not just 0.5.0)
3. Wait 5-10 minutes for Zenodo to sync
4. If still not appearing, manually upload to Zenodo

### Problem: DOI link returns 404

**Solution:**
1. Verify DOI spelling and format
2. Wait 24 hours for DOI propagation to registries
3. Check Zenodo record directly: https://zenodo.org/record/XXXXX
4. Contact Zenodo support if persistent

### Problem: Metadata not appearing in citation

**Solution:**
1. Ensure CITATION.cff is valid YAML
2. Check codemeta.json is valid JSON
3. Run `cff-converter-python` to validate:
   ```bash
   pip install cff-converter-python
   cff-convert -f cff -t bibtex CITATION.cff
   ```

---

## Additional Resources

- **Zenodo Documentation**: https://zenodo.org/help/
- **CITATION.cff Standard**: https://citation-file-format.github.io/
- **CodeMeta Standard**: https://codemeta.github.io/
- **DOI Guide**: https://www.doi.org/
- **DataCite Metadata**: https://datacite.org/
- **Force11 Software Citation**: https://www.force11.org/software-citation

---

## Next Steps

1. ✅ Create a GitHub release with proper version tag
2. ✅ Wait for Zenodo archival (automatic if integrated)
3. ✅ Update README with DOI badge
4. ✅ Share DOI in publications and documentation
5. ✅ Encourage users to cite using DOI

---

**Created**: 2024
**Last Updated**: 2024-04-13
