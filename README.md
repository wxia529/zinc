# zinc

`zinc` is a C++17 command-line tool and core library for crystal structure files. It reads common structure formats, inspects summaries, scans directories, converts files, exports optimization trajectories, and updates Quantum ESPRESSO input files.

## What It Does

- Read and auto-detect structure file formats.
- Print quick structure summaries.
- Scan directories for structure files.
- Convert between common structure formats.
- Read Quantum ESPRESSO XML, log/output, and input files.
- Update structure-related blocks in Quantum ESPRESSO inputs.

## Supported Formats

- CIF: read/write
- XYZ: read/write
- Extended XYZ trajectory: write
- VASP POSCAR/CONTCAR: read/write
- Quantum ESPRESSO XML: read
- Quantum ESPRESSO log/output: read
- Quantum ESPRESSO input: read

## Build

Requirements: CMake 3.16+ and a C++17 compiler.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Common options:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/your/install/path
cmake --install build
```

### Offline Dependencies

If the build machine cannot access the network, prepare the dependencies on a
networked machine first:

```bash
scripts/fetch_deps.sh /tmp/zinc-deps
```

Copy `/tmp/zinc-deps` to the offline machine, then configure with the generated
CMake cache file:

```bash
cmake -C /path/to/zinc-deps/zinc-deps.cmake -S . -B build
cmake --build build
```

The script downloads the exact Eigen, pugixml, CLI11, fmt, and googletest
versions used by `CMakeLists.txt`, verifies their SHA256 checksums, unpacks
them, and points `FetchContent` at those local source directories.

## Usage

Show help:

```bash
zinc --help
```

Inspect a file:

```bash
zinc info Example/yvo4.cif
zinc info Example/yvo4.cif --report-format json
zinc info Example/pwscf.in
```

Scan a directory:

```bash
zinc scan Example
zinc scan Example --element O
zinc scan Example --report-format jsonl
```

Convert formats:

```bash
zinc convert Example/yvo4.cif -o /tmp/yvo4.xyz
zinc convert Example/yvo4.cif -o /tmp/yvo4.vasp --to poscar
zinc convert Example/data-file-schema.xml -o /tmp/relax.extxyz
```

Update a Quantum ESPRESSO input file:

```bash
zinc update source.cif target.in -o updated.in
zinc update source.cif target.in --in-place
zinc update source.cif target.in --pos-unit crystal --cell-unit angstrom --fill-pseudo
```

`--from` and `--to` refer to structure formats. `--report-format` controls CLI report output. `update` supports `angstrom`, `bohr`, and `alat` for cell and position output, plus `crystal` and `crystal_sg` for positions.

## CLI Reference

### `info`

Read one structure file and print a summary.

```bash
zinc info <file> [--from <format>] [-o <output>] [--report-format text|json|jsonl]
```

### `scan`

Recursively scan a directory for structure files.

```bash
zinc scan [path] [-e|--element <element>] [--report-format text|jsonl]
```

Files without a recognized structure filename or extension are skipped.

### `convert`

Read one structure file and write it in another format.

```bash
zinc convert <file> -o <output> [--from <format>] [--to <format>] [--report-format text|json|jsonl]
```

If `--to` is omitted, the output extension decides the target format. When the target is `extxyz`, `convert` writes the full optimization trajectory if the input is a Quantum ESPRESSO output with multiple frames.

### `update`

Update the structure-related blocks in a QE input file from a source structure.

```bash
zinc update <source> <target> [-o <output>] [--in-place] [--pos-unit <unit>] [--cell-unit <unit>] [--fill-pseudo] [--report-format text|json|jsonl]
```

`--in-place` cannot be combined with `-o/--output`. If the source has no lattice, `update` tries to inherit `CELL_PARAMETERS` from the target QE input. It updates structure, atom count, and element-species information, but leaves unrelated settings such as `ecutwfc`, `occupations`, `K_POINTS`, `pseudo_dir`, and convergence controls untouched.

## Project Layout

```text
include/zinc/core/      Core data structures
include/zinc/io/        Structure driver interfaces
src/core/               Core implementations
src/io/                 File format parsers and writers
src/cli/                CLI entry point and subcommands
tests/                  GoogleTest tests
Example/                Sample input files
```

## Notes For Contributors

- Core structure type: `zinc::core::Structure`
- Format drivers implement `zinc::io::StructureDriver`
- `DriverManager` handles registration, auto-detection, reads, and writes
- To add a format, add a driver implementation and register it in `src/io/driver.cpp`
