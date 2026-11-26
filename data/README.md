# SynPEB data setup

The full dataset is not committed. Obtain SynPEB from its authors' distribution,
linked from [PPE](https://github.com/acschaefer/ppe#synthetic-plane-extraction-benchmark).
Extract it to `data/synpeb/`, preserving the original `LICENSE` and `README.md`.

Expected examples:
```
data/synpeb/train/gt/pc_01.pcd
data/synpeb/train/var1/pc_01.pcd
data/synpeb/test/gt/pc_01.pcd
data/synpeb/test/var4/pc_30.pcd
```

The archive used in this experiment is 546,383,568 bytes (about 521 MiB),
and the extracted files occupy 615,191,528 bytes. Its SHA-256 is:
```
71238e3f86e1028c79864d1ae9dc41c341317ad615b785159250262fad059a68
```
This identifies the local archive used for the report; it does not authenticate
an independently obtained download. See [synpeb.sha256](synpeb.sha256).

The downloaded archive includes a GPL version 3 LICENSE file and a very short
dataset README. The PPE software repository currently reports an MIT license;
that alone should not be treated as the dataset's redistribution terms.
Keep the original dataset attribution and license, and establish the applicable
redistribution conditions before hosting a mirror. This project links to the
authors instead of uploading a second copy.

The archive is larger than GitHub's 100 MiB regular-Git file limit. Git LFS or a
separate dataset archive can handle large distributions technically, but does
not resolve licensing questions. Keeping this input outside Git makes the code
repository smaller and preserves its provenance.

Reference: A. Schaefer et al., “A Maximum Likelihood Approach to Extract Finite
Planes from 3-D Laser Scans,” ICRA, 2019, pp. 72–78.
https://doi.org/10.1109/ICRA.2019.8794318
