# Third-party source provenance

The optional IND-CPA recomputation uses two third-party components stored in
`third_party/frodo_estimates/`.

## FrodoKEM estimates code

- Upstream repository: <https://github.com/microsoft/PQCrypto-LWEKE>
- Upstream commit: `7a4e7219d06305e16aef734213001cd8fefbcc14`
- Files used directly by this package: `cost_models.py` and the accompanying
  lattice-estimator checkout.
- License: MIT; see `third_party/frodo_estimates/LICENSE-MIT`.

The remaining small scripts and upstream README are included so that the
originating estimates directory can be inspected without another download.

## lattice-estimator

- Upstream repository: <https://github.com/malb/lattice-estimator>
- Upstream commit: `5ba00f56dd1086c3a42b98fc596c64907adb96ff`
- License stated by upstream: LGPLv3 or later.
- License texts: `LICENSE-LGPL-3.0` and `LICENSE-GPL-3.0` inside the bundled
  `lattice_estimator/` directory.

The commit is also recorded by the bundled FrodoKEM `fetch_estimator.sh`.
All third-party source is kept separate from the LeechKEM reproduction code.
