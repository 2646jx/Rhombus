# Rhombus: Fast Homomorphic Matrix-Vector Multiplication for Secure Two-Party Inference

## Introduction

This repository provides the implementation of the paper [Rhombus](https://eprint.iacr.org/2024/1611.pdf) (except the network communication).
In [Rhombus](https://eprint.iacr.org/2024/1611.pdf), the authors propose the matrix-vector multiplication (MVM) and matrix-matrix multiplication (MatMul).
For MVM, there are row-major (based on PackRLWEs) and column-major (based on Expand) approaches, we now only implement the row-major approach, the column-major
approach might be updated later. For the MatMul algorithm, [Rhombus](https://eprint.iacr.org/2024/1611.pdf) proposed the PackRLWEs, Expand based approaches,
and the V1, V2 splitting method for matrix partitioning. In summary, there are 4 algorithms: PackRLWEs+V1, PackRLWEs+V2, Expand+V1, Expand+V2.
This project implements three algorithms: PackRLWEs+V1, PackRLWEs+V2, and Expand+V2.

## Building

`scripts/`, `src/`, `test/` are the main directories of this project.

- `scripts/`: there are two scripts `build-deps.sh` and `build.sh`, which are used to build the dependency [SEAL](https://github.com/microsoft/SEAL) and to build this project.
- `src/`: the implementation of [Rhombus](https://eprint.iacr.org/2024/1611.pdf)
- `test/`: the module tests

To build this project, execute the following commands:

```PowerShell
bash scripts/build-deps.sh
bash scripts/build.sh
```

then the library and tests will be built.

## Test

Once built successfully, run the following two commands to test the matrix-vector multiplication and matrix-matrix multiplication, respectively:

```PowerShell
./build/bin/matvec
```

```PowerShell
./build/bin/matmul
```

## Options

[SEAL](https://github.com/microsoft/SEAL) library has two useful options to improve the performance of homomorphic computation.
To config the two options, you can open the scripts file `scripts/build-deps.sh`:

- `USE_HEXL`: set to `ON` to enable the hexl acceleration provided that your platform supports the Intel AVX512 instructions.
- `USE_ZSTD`: set to `ON` to enable the zstd compression, which can compress the ciphertexts transmitted.

After resetting the options, execute the commands again:

```PowerShell
bash scripts/build-deps.sh
bash scripts/build.sh
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
