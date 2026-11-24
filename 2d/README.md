# 二维求解器

## 编译与运行

在 Linux 或 WSL 的 Bash 中，从 `2d` 目录运行。

```bash
make
./ugks2d cases/poiseuille_fast.namelist
```

Windows PowerShell 使用 GNU Fortran 编译并启动同一算例。

```powershell
./build.ps1
./ugks2d.exe cases/poiseuille_fast.namelist
```

## 选择算例

| 要运行的算例 | 输入文件或命令 |
|---|---|
| Poiseuille 流动 | `cases/poiseuille_fast.namelist` |
| Rayleigh–Taylor 不稳定性 | `cases/rt_60_fast.namelist` |
| 单组分或五组分轴对称流动 | `./ugks2d --axisymmetric CASE.nml`，见[轴对称使用说明](../docs/axisymmetric-quickstart.md) |

不带输入参数运行时，程序读取当前目录的 `input.namelist`。二维笛卡尔算例输出 `macro*.dat`。轴对称输出由输入中的物种配置和平均时段决定。
