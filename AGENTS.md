# Repository Guidelines

## Project Structure & Module Organization

`zinc` is a C++17 command-line tool and core library for crystal structure files.

- `include/zinc/core/`: public core types such as `Lattice`, `Structure`, and constants.
- `include/zinc/io/`: public structure driver interfaces and format driver declarations.
- `src/core/`: core implementations.
- `src/io/`: format parsers/writers for CIF, XYZ, POSCAR, Quantum ESPRESSO XML/log/input, and driver registration.
- `src/cli/`: CLI entry point and subcommands: `info`, `scan`, `convert`, and `update`.
- `tests/`: GoogleTest unit and CLI tests.
- `Example/`: sample structure, QE input/output, and trajectory files used by tests and manual checks.

## Build, Test, and Development Commands

```bash
cmake -S . -B build
```
Configure the project and fetch dependencies with CMake `FetchContent`.

```bash
cmake --build build
```
Build the `zinc` CLI and test binary.

```bash
ctest --test-dir build --output-on-failure
```
Run the full GoogleTest suite.

```bash
./build/zinc info Example/test.cif
./build/zinc convert Example/data-file-schema.xml -o /tmp/relax.extxyz
```
Run local CLI smoke checks.

## Coding Style & Naming Conventions

Use C++17 and follow the existing style: two-space indentation, braces on the same line for functions and control blocks, and concise helper functions in anonymous namespaces where appropriate. Public classes use `PascalCase`; functions and variables use `snake_case`. Keep comments sparse and focused on non-obvious parsing or unit-conversion logic.

When adding a new format, implement `zinc::io::StructureDriver`, add the public header under `include/zinc/io/`, the implementation under `src/io/`, register it in `src/io/driver.cpp`, and update `CMakeLists.txt`.

## Testing Guidelines

Tests use GoogleTest. Add focused tests near the behavior being changed:

- core behavior in `tests/test_core.cpp`
- parser/writer behavior in `tests/test_io.cpp`
- command-line behavior in `tests/test_cli.cpp`

Prefer fixture files in `Example/` for realistic parser coverage. For generated temporary files, write to `std::filesystem::temp_directory_path()` and remove them after the test.

## Commit & Pull Request Guidelines

This checkout does not include Git history, so no repository-specific commit convention is available. Use short imperative commit messages, for example `Add QE input parser` or `Fix extxyz frame count`.

Pull requests should include a clear summary, affected commands or formats, test results, and any compatibility notes. For parser changes, mention representative input files and edge cases covered.

## Agent-Specific Instructions

Do not rewrite unrelated files or reformat broad areas. Keep parser changes conservative, preserve existing CLI behavior where possible, and update README/tests when user-facing behavior changes.
