# GEO download troubleshooting

If `GEOquery::getGEO()` fails with an error such as:

```r
Failed to download data/raw/GEO/GSE25097_series_matrix.txt.gz
```

the most common cause is that the GEO series does not expose a simple
`<GSE>_series_matrix.txt.gz` file or the download times out. The revised
`R/_shared.R` implements a fallback that:

1. lists the official NCBI GEO matrix directory;
2. detects the actual `*series_matrix.txt.gz` file name(s);
3. downloads them directly with retries and a longer timeout;
4. parses the local series matrix with `GEOquery::getGEO(filename=...)`;
5. fetches the GPL annotation separately if the fallback matrix lacks gene symbols.

Run again from the repository root:

```bash
Rscript scripts/02_download_and_curate_geo.R
```

If a partially downloaded zero-byte file exists, delete it first:

```bash
rm -f data/raw/GEO/GSE25097*
```
