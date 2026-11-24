# Hermite-UGKS--Axisymmetric

## 编译

在 Linux 或 WSL 中安装 GNU Fortran 和 GNU Make，然后运行。

```bash
git clone https://github.com/zbyandmoon/Hermite-UGKS-axisym.git
cd Hermite-UGKS-axisym
make -C 2d
```

Windows 使用 PowerShell 和 GNU Fortran，运行 `./2d/build.ps1` 生成 `2d/ugks2d.exe`。

## 运行

单组分和五组分都使用二维可执行文件的 `--axisymmetric` 入口。单组分将 `species_file` 留空，五组分指定包含五种物性的文件。编译、完整运行命令、参数和输出文件见[轴对称使用说明](docs/axisymmetric-quickstart.md)。

其他求解入口见[一维用法](1d/README.md)和[二维用法](2d/README.md)。
