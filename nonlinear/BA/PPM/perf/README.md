# Microarchitecture Performance

## Model Table
| Compute Node | Microarchitecture | Cores/Node | PBS Code |
|--------------|-------------------|------------|----------|
| Aitken       | Cascade Lake      | 40         | cas_ait  |
| Electra      | Skylake           | 40         | sky_ele  |
| Pleiades     | Broadwell         | 28         | bro      |
|              | Haswell           | 24         | has      |
|              | Ivy Bridge        | 20         | ivy      |
|              | Sandy Bridge      | 16         | san      |

_From [Basic Tasks: Submitting Interactive Jobs (HECC KB)](https://www.nas.nasa.gov/hecc/support/kb/basic-tasks_264.html#Submitting%20Interactive%20Jobs)_

## [Athena++ Configuration](https://github.com/PrincetonUniversity/athena/wiki/Configuring#options)
To specify a target architecture, add flags (overriding incompatible, earlier flags generated by the configure script) with the `--cflag` option.
Note this option should be set with an `=` sign to avoid the script trying to interpret any dashes, e.g., `--cflag="-axCORE-AVX512,CORE-AVX2 -xAVX"`.

## HECC KB Resources
- [**Recommended Compiler Options**](https://www.nas.nasa.gov/hecc/support/kb/recommended-compiler-options_99.html)
- Preparing to run _only_ on
  - [Pleiades Sandy Bridge Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-pleiades-sandy-bridge-nodes_322.html)
  - [Pleiades Ivy Bridge Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-pleiades-ivy-bridge-nodes_446.html)
  - [Pleiades Haswell Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-pleiades-haswell-nodes_491.html)
  - [Pleiades Broadwell Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-pleiades-broadwell-nodes_530.html)
  - [Electra Skylake Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-electra-skylake-nodes_551.html)
  - [Aitken Cascade Lake Nodes](https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-aitken-cascade-lake-nodes_597.html)
- Configuration Details
  - [Pleiades](https://www.nas.nasa.gov/hecc/support/kb/pleiades-configuration-details_77.html)
  - [Electra](https://www.nas.nasa.gov/hecc/support/kb/electra-configuration-details_537.html)
  - [Aitken](https://www.nas.nasa.gov/hecc/support/kb/aitken-configuration-details_580.html)
