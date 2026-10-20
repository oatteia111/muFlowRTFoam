<#
.SYNOPSIS
Runs muFlowRT with the DLL paths for a local UCRT64 build or an extracted runtime ZIP.
.EXAMPLE
& .\scripts\run-muflowrt.ps1 -SolverArguments '-help'
.EXAMPLE
& .\scripts\run-muflowrt.ps1 -CaseDirectory 'C:\models\my-case'
#>
param(
    [string]$CaseDirectory = (Get-Location).Path,
    [string]$OpenFOAMRoot = 'C:\OF13\OpenFOAM-13-Windows',
    [string]$PhreeqcRMRoot = 'C:\phreeqcrm\phreeqcrm-3.8.6-17100',
    [string]$Msys2Root = 'C:\msys64',
    [string]$ExecutablePath = 'C:\OF13User\of-user\platforms\mingw_w64GccDPInt32Opt\bin\muFlowRT.exe',
    [string[]]$SolverArguments = @()
)

$ErrorActionPreference = 'Stop'
# Preserve the native exit code even if the caller enables this PowerShell 7 option.
$PSNativeCommandUseErrorActionPreference = $false
$previousEnvironment = @{}
$exitCode = 1

try {
    $wmOptions = 'mingw_w64GccDPInt32Opt'
    $bundleBin = Join-Path $PSScriptRoot 'bin'
    $bundleFoam = Join-Path $PSScriptRoot 'OpenFOAM-13'
    if ((Test-Path -LiteralPath $bundleBin -PathType Container) -or
        (Test-Path -LiteralPath $bundleFoam -PathType Container)) {
        $OpenFOAMRoot = $bundleFoam
        $ExecutablePath = Join-Path $bundleBin 'muFlowRT.exe'
        $foamLibBin = $bundleBin
        $pstreamBin = $bundleBin
        $phreeqcBin = $bundleBin
        $runtimeBin = $bundleBin
        $userLibBin = $bundleBin
    }
    else {
        # Resolve overrides before changing to the case directory.
        $OpenFOAMRoot = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OpenFOAMRoot)
        $PhreeqcRMRoot = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($PhreeqcRMRoot)
        $Msys2Root = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Msys2Root)
        $ExecutablePath = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($ExecutablePath)
        $foamLibBin = Join-Path $OpenFOAMRoot "platforms\$wmOptions\lib"
        $pstreamBin = Join-Path $foamLibBin 'dummy'
        $phreeqcBin = Join-Path $PhreeqcRMRoot 'bin'
        $runtimeBin = Join-Path $Msys2Root 'ucrt64\bin'
        $userLibBin = Join-Path (Split-Path (Split-Path $ExecutablePath -Parent) -Parent) 'lib'
    }

    $requiredFiles = @(
        $ExecutablePath
        (Join-Path $OpenFOAMRoot 'etc\controlDict')
        (Join-Path $pstreamBin 'libPstream.dll')
        (Join-Path $phreeqcBin 'libPhreeqcRM.dll')
    )
    foreach ($library in 'OpenFOAM', 'finiteVolume', 'meshTools', 'fileFormats', 'surfMesh', 'triSurface') {
        $requiredFiles += Join-Path $foamLibBin "lib$library.dll"
    }
    foreach ($dll in 'libgcc_s_seh-1.dll', 'libstdc++-6.dll', 'libwinpthread-1.dll', 'libgomp-1.dll', 'zlib1.dll') {
        $requiredFiles += Join-Path $runtimeBin $dll
    }
    $missing = @($requiredFiles | Where-Object { -not (Test-Path -LiteralPath $_ -PathType Leaf) })
    if ($missing.Count) {
        throw "Missing runtime files. Check the local path parameters or extract the complete ZIP:`n$($missing -join "`n")"
    }

    $ExecutablePath = (Resolve-Path -LiteralPath $ExecutablePath).ProviderPath
    $OpenFOAMRoot = (Resolve-Path -LiteralPath $OpenFOAMRoot).ProviderPath
    $CaseDirectory = (Resolve-Path -LiteralPath $CaseDirectory).ProviderPath
    if (($SolverArguments -notcontains '-help') -and
        -not (Test-Path -LiteralPath (Join-Path $CaseDirectory 'system\controlDict') -PathType Leaf)) {
        throw "No system/controlDict in $CaseDirectory. Use -CaseDirectory with a prepared muFlowRT case, or -SolverArguments '-help' to check startup."
    }

    $dllDirectories = @($foamLibBin, $pstreamBin, $userLibBin, $phreeqcBin, $runtimeBin)
    $environment = @{
        WM_PROJECT = 'OpenFOAM'
        WM_PROJECT_VERSION = '13'
        WM_PROJECT_DIR = $OpenFOAMRoot.Replace('\', '/')
        WM_PROJECT_INST_DIR = (Split-Path $OpenFOAMRoot -Parent).Replace('\', '/')
        WM_OPTIONS = $wmOptions
        WM_MPLIB = 'Dummy'
        FOAM_MPI = 'dummy'
        FOAM_ETC = (Join-Path $OpenFOAMRoot 'etc').Replace('\', '/')
        FOAM_LIBBIN = $foamLibBin.Replace('\', '/')
        FOAM_USER_LIBBIN = $userLibBin.Replace('\', '/')
        FOAM_APPBIN = (Split-Path $ExecutablePath -Parent).Replace('\', '/')
        FOAM_USER_APPBIN = (Split-Path $ExecutablePath -Parent).Replace('\', '/')
        PATH = (@($dllDirectories | Select-Object -Unique) + "$env:SystemRoot\System32" + $env:SystemRoot) -join ';'
    }
    foreach ($name in $environment.Keys) {
        $previousEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
        [Environment]::SetEnvironmentVariable($name, $environment[$name], 'Process')
    }

    Push-Location -LiteralPath $CaseDirectory
    try {
        Write-Host "muFlowRT: $ExecutablePath"
        Write-Host "Case directory: $CaseDirectory"
        & $ExecutablePath @SolverArguments
        $exitCode = $LASTEXITCODE
        if ($exitCode -ne 0) {
            Write-Warning "muFlowRT exited with code $exitCode."
        }
    }
    finally {
        Pop-Location
    }
}
catch {
    Write-Error $_ -ErrorAction Continue
}
finally {
    foreach ($name in $previousEnvironment.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previousEnvironment[$name], 'Process')
    }
}
exit $exitCode
