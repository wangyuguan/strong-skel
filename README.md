# strong-skel

## HIRAX compatibility changes in this fork

This fork starts from upstream commit
`e64dfdcb68186b5545910d5e7002c5765d107145`. The four MATLAB runtime files
`srskelf_asym.m`, `shypoct.m`, `srskelf_sv_nn.m`, and `srskelf_mv_nn.m`
were synchronized from the HIRAX server dependency on 2026-10-11.
The upstream authorship and GPL-2.0-or-later license are retained.

- `srskelf_asym(..., [], opts)` accepts an empty proxy callback. At levels
  greater than 2, the ID samples **every active non-neighbor DOF**, in both
  matrix directions for asymmetric input, including accumulated Schur updates.
  It uses FLAM's deterministic `id`; it does not use randomized sampling.
- Near-field interactions remain in the elimination. Level 2 uses the
  upstream weak skeletonization step, and the root uses a dense local factor.
- The non-transposed asymmetric multiply and solve routines are also available
  in the repository root, so HIRAX can use `addpath(strong_skel_root)`.
  The original `mv/` and `sv/` files are retained.
- A nonempty proxy callback keeps the upstream behavior. Other factorization
  variants have not been changed. Existing `opts.verb` progress output is retained.

No full original matrix is assembled by this change, but each box's sampled
far-field block is formed explicitly. This mode can therefore be more expensive
than a proxy-based factorization.

For a small regression test, add FLAM to the MATLAB path, add this repository's
`test` directory, and run `test_empty_proxy`. It checks all-far sampling on a
four-level tree, actual elimination, complex asymmetric multiplication and
solves, and preservation of the random-number state.


A research implementation of RS-S and RS-WS from Minden et al., "A recursive skeletonization factorization based on strong admissibility" ([Journal version](http://dx.doi.org/10.1137/16M1095949), [arXiv](https://arxiv.org/abs/1609.08130))

We emphasize that this code is provided under the GNU GPL as an as-is research reference implementation.  It is not an actively maintained library.  

Contributors:
Victor Minden, Ken L. Ho, Anil Damle, Lexing Ying

## Details
Built on the [FLAM library](https://github.com/klho/FLAM/) by Ken L. Ho, this research code implements the factorizations described in "A recursive skeletonization factorization based on strong admissibility" by Minden et al.

The main functions provided are as follows:
- `srskelf.m`: the strong recursive skeletonization factorization (aka "RS-S") for symmetric matrices.  This function uses Cholesky factorizations of diagonal subblocks that can be assumed positive-definite in exact arithmetic.
- `srskelf_asym.m`: RS-S for asymmetric matrices.  This function uses LU factorizations where `srskelf.m` uses Cholesky.
- `srskelf_hybrid.m`: the hybrid recursive skeletonization factorization (aka "RS-WS") for symmetric matrices.  This function is like `srskelf.m`, but additionally interlaces levels of standard ("weak") skeletonization.
- `srskelf_asymhybrid.m`: RS-WS for asymmetric matrices.  This function is like `srskelf_asym.m`, but additionally interlaces levels of standard ("weak") skeletonization.

Tests are available as MATLAB functions in the directory "test" and may be run with no arguments, assuming that FLAM and the strong-skel directories are in your MATLAB path.  The tests provided are as follows:
- `ie_square.m`: The setup of Example 1 as described in Minden et al.
- `ie_cube.m`: The setup of Example 2 as described in Minden et al.
- `ie_sphere.m`: The setup of Example 3 as described in Minden et al.