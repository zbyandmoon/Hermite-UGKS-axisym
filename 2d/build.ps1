$ErrorActionPreference = 'Stop'

$sourceDir = Join-Path $PSScriptRoot 'src'
$buildDir = Join-Path $PSScriptRoot 'build'
$executable = Join-Path $PSScriptRoot 'ugks2d.exe'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Push-Location -LiteralPath $buildDir

try {
    $flags = @('-O2','-ffree-line-length-none','-fopenmp')
    $modules = @(
        'config','data2d','grid','gastheory2d','gravity','hermite_symmetric2d','velocity',
        'boundary','flux','initialization','solver','output2d',
        'axisymmetric_kinetic','axisymmetric_mesh','axisymmetric_internal_modes','axisymmetric_mixture_collision',
        'axisymmetric_heat_flux','axisymmetric_plume','axisymmetric_mixture','main'
    )

    foreach ($name in $modules) {
        & gfortran @flags '-J.' '-I.' '-c' (Join-Path $sourceDir "$name.f90") '-o' "$name.o"
        if ($LASTEXITCODE -ne 0) {
            throw "Compilation failed for $name.f90"
        }
    }

    $objects = $modules | ForEach-Object { "$_.o" }
    # Link in the build directory so GNU ld does not need a Unicode output path.
    & gfortran '-fopenmp' @objects '-o' 'ugks2d.exe'
    if ($LASTEXITCODE -ne 0) {
        throw 'Link failed'
    }
    Move-Item -Force -LiteralPath '.\ugks2d.exe' -Destination $executable
    Write-Host "Built $executable"
}
finally {
    Pop-Location
}
