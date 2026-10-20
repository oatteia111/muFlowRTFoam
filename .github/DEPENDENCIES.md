# Windows muFlowRT builds

One manual workflow, **Build muFlowRT Windows**, builds PhreeqcRM and the seven
OpenFOAM libraries concurrently, then builds and publishes muFlowRT. Its unsigned
release contains a Windows runtime ZIP and both reusable SDKs. Anchor QEA signing
and JFrog publication remain separate work.

## Local builds

Use an MSYS2 **UCRT64** terminal. From this repository:

```sh
set -o pipefail
bash scripts/build-phreeqcrm.sh 2>&1 | tee phreeqcrm-build.log
export WM_NCOMPPROCS=4
bash scripts/build-openfoam.sh 2>&1 | tee openfoam-limited-build.log

# Once both dependencies are available:
bash scripts/build-muflowrt.sh 2>&1 | tee muflowrt-build.log
```

The solver script locates the repository relative to itself, so it can also be
invoked by absolute path from another directory. It consumes existing SDKs,
selects dummy Pstream, applies Olivier's `LIB_DIR`/`LIB_HEADER_DIRS` exports, and
runs `wmake` directly in `solvers`. It prints the executable location, normally
`C:\OF13User\of-user\platforms\mingw_w64GccDPInt32Opt\bin\muFlowRT.exe`.

| Solver environment variable | Default |
|---|---|
| `OF13_ROOT` | `/c/OF13` |
| `OF13_CLONE` | `$OF13_ROOT/OpenFOAM-13-Windows` |
| `PHREEQCRM_ROOT` | `/c/phreeqcrm/phreeqcrm-3.8.6-17100` |
| `OF13_USER` | Upstream default `/c/OF13User/of-user`; preserve caller override |
| `WM_NCOMPPROCS` | Preserve caller value; otherwise upstream uses `nproc` |

Use short MSYS2 paths without spaces for SDKs and source checkouts, as required
by wmake's include/link recipes. Overrides apply to the solver; dependency
builders retain their fixed installation locations.

An interrupted build can reuse its objects and libraries when sources and the
toolchain are unchanged. No deletion is needed before the first solver build.
The upstream dependency scanner can miss header dependencies; header edits may
require a clean solver rebuild until that upstream issue is fixed.

### Run a locally built solver

The build script's environment does not persist in the parent terminal. Running
the EXE directly in PowerShell therefore cannot find the OpenFOAM/PhreeqcRM DLLs.
Use `scripts/run-muflowrt.ps1` from an **existing PowerShell terminal**:

```powershell
$launcher = 'C:\Users\mgreen\git\muFlowRTFoam\scripts\run-muflowrt.ps1'
& $launcher -SolverArguments '-help'
& $launcher -CaseDirectory 'C:\models\my-case'
```

The help command checks startup without running a simulation. Replace the case
path with a prepared muFlowRT case. The default case directory is your current
directory. Use `-SolverArguments @('-option', 'value')` for solver arguments.
The launcher temporarily sets the OpenFOAM environment and DLL paths, then
restores the caller's environment/directory and returns the solver's exit code.
For custom local installations, override `-OpenFOAMRoot`, `-PhreeqcRMRoot`,
`-Msys2Root`, or `-ExecutablePath` using Windows paths.

If execution policy blocks the unsigned script, allow it in this terminal only:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

Run the launcher from an existing PowerShell terminal so the output remains
visible after the solver exits.

Alternatively, before running the executable from an MSYS2 UCRT64 terminal:

```sh
export OF13_ROOT="${OF13_ROOT:-/c/OF13}"
export OF13_CLONE="${OF13_CLONE:-$OF13_ROOT/OpenFOAM-13-Windows}"
export PHREEQCRM_ROOT="${PHREEQCRM_ROOT:-/c/phreeqcrm/phreeqcrm-3.8.6-17100}"
export WM_MPLIB=Dummy FOAM_MPI=dummy
source "$OF13_CLONE/scripts/windows/env.sh"
export LIB_DIR="$WM_PROJECT_DIR/src"
export LIB_HEADER_DIRS="$WM_PROJECT_DIR/src"
export PATH="$PHREEQCRM_ROOT/bin:$PATH"

cd /c/path/to/your/case
"$FOAM_USER_APPBIN/muFlowRT.exe"
```

The solver loads `phreeqc.dat` from its working directory. Put the intended
database in the case directory. The SDK supplies a reference copy under
`$PHREEQCRM_ROOT/database/`; do not overwrite a case's chosen database without
checking it.

## Toolchain and dependency recipes

Each CI job provisions its own runner through `.github/actions/setup-msys2`.
The shared base installs UCRT64 GCC, libgomp, zlib, make, Git, curl, tar, and zstd.
The action accepts an optional `extra-packages` input:

| Job | Additional packages |
|---|---|
| PhreeqcRM | UCRT64 CMake and Ninja |
| OpenFOAM | flex and bison |
| muFlowRT | zip; wmake comes from the OpenFOAM SDK |

MSYS2 updates and package caching stay enabled. Each job records its full package
list and a smaller compiler/runtime record. The final job requires both SDKs to
match its GCC version, architecture, and runtime package versions. Different
build-tool packages are expected. If rolling updates cause a mismatch, rerun all
jobs together; do not bypass the compatibility check.

### PhreeqcRM

`scripts/build-phreeqcrm.sh` downloads **3.8.6-17100** from GitHub release `v3.8.6`.
It uses native UCRT64 CMake/Ninja, Release, shared libraries, and OpenMP. MPI,
Python, Fortran testing, YAML, and zlib are disabled. OpenMP requires
`mingw-w64-ucrt-x86_64-libgomp`, installed by the shared action.

Downloads/sources go under `/c/phreeqcrm/source`, compilation under
`/c/phreeqcrm/build`, and installation under `/c/phreeqcrm/phreeqcrm-3.8.6-17100`.
The SDK preserves headers, `irm_dll_export.h`, import library, DLL, CMake metadata,
databases, documentation, notices, source archive, and CMake build options.
The source SHA-256 is:

```text
d2e532a065b17ab07104c7f445c823aadafab10b263b613b1fec7e5fb50b046a
```

### OpenFOAM

`scripts/build-openfoam.sh` preserves the hardcoded OpenFOAM and ThirdParty
commits captured from their upstream default branches on 2026-10-19. Change these
pins in the script to change sources.

Sibling checkouts live under `C:\OF13`. The script builds wmake tools, Windows
support objects, the core/Pstream circular-link bootstrap, and these libraries:
**OpenFOAM, dummy Pstream, fileFormats, surfMesh, triSurface, meshTools,
finiteVolume**. ThirdParty is retained as source; its libraries and OpenFOAM
applications such as checkMesh are not built.

The SDK includes both source trees, generated headers, materialized lnInclude
contents, wmake/support tools, DLLs/import libraries, configuration, and licenses.
Git metadata, object/dependency files, and build logs are excluded. SDK consumers
need a matching UCRT64 toolchain; the runtime ZIP supplies the DLLs needed to run
the solver without that toolchain.

## Workflow and release assets

The workflow definition must exist on the default branch, `main`, before manual
dispatch is available. Select `cicd` to run that branch's implementation:

```sh
gh workflow run build-muflowrt.yml --ref cicd
```

Dependency SDKs/checksums and diagnostics have **14-day** Actions artifact
retention. Job outputs carry artifact names so the final job can reuse successful
dependency builds when rerun. It verifies checksums, restores SDKs, builds the
solver, and uploads release assets before publication. Only the final job has
`contents: write`.

The prerelease tag is `muflowrt-ucrt64-x86_64-<run-id>-<attempt>`, targeting the
workflow's exact repository commit. It is never marked Latest. Assets are:

- `muflowrt-ucrt64-x86_64.zip`: executable, seven OpenFOAM DLLs, PhreeqcRM DLL,
  compiler runtimes, OpenFOAM configuration, reference database, launcher,
  notices, build metadata, and internal checksums.
- `phreeqcrm-3.8.6-17100-ucrt64-x86_64.tar.zst`: installed PhreeqcRM SDK.
- `openfoam-13-ucrt64-x86_64.tar.zst`: OpenFOAM and ThirdParty SDK trees.
- `BUILD_INFO.txt`, `TOOLCHAIN_ABI.txt`, `SDK-README.md`, and `SHA256SUMS`.

Packaging checks imports of every included executable/DLL against packaged DLLs
and known Windows system libraries. The workflow does not run smoke tests,
upstream test suites, or numerical validation.

### Run the runtime ZIP

Extract the whole ZIP and, from a prepared case directory, run:

```powershell
& 'C:\path\to\muflowrt-ucrt64-x86_64\run-muflowrt.ps1'
```

The PowerShell launcher detects the extracted bundle layout. It also accepts
`-CaseDirectory 'C:\models\my-case'` and `-SolverArguments '-help'`, as above.
The launcher forwards arguments/exit status and uses bundled DLLs plus Windows
system directories in PATH. The case still needs its intended `phreeqc.dat`.
Normal serial cases need no MSYS2 installation. Runtime
compilation (`#codeStream`/`#calc`), extra plugins, MPI, and other OpenFOAM utilities
require additional tools outside this runtime bundle.

### Extract SDKs for future builds

Download all assets from an explicit release tag into an empty directory.
In an MSYS2 UCRT64 shell:

```sh
sha256sum --check --strict SHA256SUMS
mkdir -p /c/phreeqcrm
zstd -dc phreeqcrm-3.8.6-17100-ucrt64-x86_64.tar.zst | tar -xf - -C /c/phreeqcrm
zstd -dc openfoam-13-ucrt64-x86_64.tar.zst | tar -xf - -C /c
```

Archives exceeding **1,900 MiB** are split into ordered parts. Verify their
checksums, then reassemble during extraction:

```sh
cat openfoam-13-ucrt64-x86_64.tar.zst.part-* | zstd -dc | tar -xf - -C /c
# Use the same pattern for split PhreeqcRM, extracting to /c/phreeqcrm.
# Reassemble a split runtime ZIP before opening it:
cat muflowrt-ucrt64-x86_64.zip.part-* > muflowrt-ucrt64-x86_64.zip
```

Source OpenFOAM's environment and apply the header exports in each downstream
build step that needs them, as shown in the local-run example.

### Publication retries

A draft is published only after all asset uploads complete. Identical existing
assets are reused; conflicting replacements are rejected.

Rerunning the failed final job reuses dependency SDKs, rebuilds the solver, and
uses the new attempt number. To retry only publication, download the retained
final release-assets artifact and run this in Git Bash with an authenticated gh:

```sh
export GITHUB_REPOSITORY=AnchorQEA/muFlowRTFoam
export GITHUB_SHA=WORKFLOW_COMMIT_FROM_BUILD_INFO
bash scripts/publish-dependency.sh ORIGINAL_RELEASE_TAG /c/path/to/release-assets
```

Packaging helpers use CI staging variables; the three build scripts are the local
entry points.

## References

- [OpenFOAM Windows recipe](https://github.com/fralarco/OpenFOAM-13-Windows/blob/master/BUILD_WINDOWS.md)
- [PhreeqcRM release](https://github.com/phreeqc-dev/phreeqcrm/releases/tag/v3.8.6)
- [MSYS2 setup action](https://github.com/msys2/setup-msys2)
- [GitHub workflow artifacts](https://docs.github.com/en/actions/tutorials/store-and-share-data)
