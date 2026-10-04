# Full KEM runtime results

Times are milliseconds per operation and shown as median (p95).

| Scheme | KeyGen | Encaps | Decaps valid | Decaps invalid | Full cycle |
|---|---:|---:|---:|---:|---:|
| LeechKEM-1 | 33.759 (41.128) | 13.878 (14.314) | 14.280 (17.701) | 14.205 (14.749) | 61.409 (63.700) |
| FrodoKEM-640-SHAKE | 5.317 (6.054) | 9.188 (9.509) | 9.017 (9.284) | 9.013 (9.229) | 23.454 (24.323) |
| LeechKEM-2 | 54.569 (60.119) | 23.670 (25.562) | 24.170 (25.140) | 24.232 (25.307) | 102.495 (109.023) |
| FrodoKEM-976-SHAKE | 11.661 (12.540) | 17.474 (18.357) | 17.313 (18.027) | 17.423 (18.402) | 46.099 (47.063) |
| LeechKEM-3 | 114.375 (134.015) | 47.738 (50.528) | 48.449 (51.571) | 48.436 (50.309) | 209.628 (228.977) |
| FrodoKEM-1344-SHAKE | 21.209 (22.061) | 36.977 (38.389) | 36.642 (37.693) | 36.486 (37.052) | 95.693 (102.150) |

## Paired median ratios (LeechKEM / FrodoKEM)

| Pair | KeyGen | Encaps | Decaps valid | Decaps invalid | Full cycle |
|---|---:|---:|---:|---:|---:|
| LeechKEM-1 / FrodoKEM-640-SHAKE | 6.35x | 1.51x | 1.58x | 1.58x | 2.62x |
| LeechKEM-2 / FrodoKEM-976-SHAKE | 4.68x | 1.35x | 1.40x | 1.39x | 2.22x |
| LeechKEM-3 / FrodoKEM-1344-SHAKE | 5.39x | 1.29x | 1.32x | 1.33x | 2.19x |
