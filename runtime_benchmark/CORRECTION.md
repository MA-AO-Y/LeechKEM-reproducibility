# g_tau indexing correction

The preceding benchmark package grouped the `tau-1` residual bits of each
coordinate contiguously.  The manuscript instead orders these bits by bit
plane: coordinate `j` uses `b_{g,24r+j}` for `0 <= r <= tau-2`.

This release changes both directions in `src/leechkem_ref.c`:

- `lattice_encode` reads each residual digit using `24*r+j`;
- `lattice_decode` writes each recovered digit using `24*r+j`.

The self-test independently sets every `b_{g,24r+j}` basis bit and checks that
the encoded vector has value `2^(r+1)*beta` in coordinate `j` and zero in all
other coordinates.  The existing exhaustive message-basis round trips and the
random round trips then check the inverse mapping.

Because the change affects the timed implementation, paper measurements must
be regenerated with this release on the reporting machine.
