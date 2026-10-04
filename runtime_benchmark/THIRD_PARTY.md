# Third-party code

This package contains source snapshots needed for a reproducible offline run.

- `third_party/avanpo_leech`: Vardy--Be'ery Leech-lattice decoder from
  `avanpo/leech-decoding`, commit
  `fefc3e3fe14819195297cd6f7c81348fd9965166`.  License: GNU GPL v3; see the
  adjacent `LICENSE` file.  The decoding algorithm is unchanged.  One entry
  in the separate `leech_utils.c` MOG generator table is corrected: row 19 is
  restored to the standard row printed in Figure 3.1 of the author's thesis.
- `third_party/PQCrypto-LWEKE`: Microsoft Research FrodoKEM reference code and
  its SHAKE implementation from `microsoft/PQCrypto-LWEKE`, commit
  `e1edeb3af1fae0d5727683bd2f5465280ec2437a`.  License: MIT; see the adjacent
  `LICENSE` file.

The top-level package is distributed under GNU GPL v3 because it links the
GPL-licensed Leech decoder.
