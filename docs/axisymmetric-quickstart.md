# 轴对称 Hermite-UGKS 快速入门

同一个二维可执行文件通过 `--axisymmetric` 进入轴对称模式。`species_file` 为空时运行单组分。指定含五行物性的文件时运行五组分。多组分入口也接受一行等效物性，可用于在相同数值选项下做单组分对照。

## 编译

在 Linux 或 WSL 的 Bash 中，从仓库根目录运行。需要 GNU Fortran、GNU Make 和 OpenMP。

```bash
make -C 2d
```

在装有 `gfortran` 的 Windows PowerShell 中，也可运行 `./2d/build.ps1`，得到 `2d/ugks2d.exe`。

## 运行两种示例

以下是封闭圆柱腔内的短算例，初态为 300 K 静止均匀气体。五组分采用示例物性和二阶 Minmod。两个算例分别写入独立目录。

```bash
repo=$PWD
export OMP_NUM_THREADS=2
mkdir -p runs
out=$(mktemp -d "$repo/runs/axisymmetric.XXXXXX")
mkdir "$out/single" "$out/five"

(
    cd "$out/single"
    "$repo/2d/ugks2d" --axisymmetric "$repo/2d/cases/axisymmetric/single-species.nml" > run.log 2>&1
)

cp "$repo/2d/cases/axisymmetric/five-species.dat" "$out/five/"
(
    cd "$out/five"
    "$repo/2d/ugks2d" --axisymmetric "$repo/2d/cases/axisymmetric/five-species.nml" > run.log 2>&1
)
printf '%s\n' "$out"
```

两个示例均推进至 35 µs，从 10 µs 起统计。完成时 `run.log` 给出实际平均时长。单组分平均场为 `axisymmetric_average.dat`，五组分为 `axisymmetric_mixture_average.dat`。后者同时输出各组分密度、实际轴向与径向网格面，以及逐组分质量收支。

## 配置自己的算例

- 在 `&axisym` 中设置物理域、网格、初始状态、壁温、`end_time` 与 `average_start`。长度用 m，时间用 s，温度用 K。
- 喷口流动采用 `uniform_test=.false.`，设置机体位置 `x_front`、`x_aft`、`r_body` 和径向入口文件 `inlet_file`。入口表首行为表头，后续各行依次给出半径下界、半径上界、密度、温度、轴向速度和径向速度。文件路径相对于运行目录解析。
- 五组分表的列为物种名、摩尔质量（g/mol）、质量分数、VHS 参考直径（m）、黏度指数、转动自由度和转动交换概率。质量分数按所用工况设置。`vhs_reference_temperature` 应与物种直径的参考温度一致。
- 多组分路径用 `spatial_order=2` 启用二阶重构，`slope_limiter=1` 为 Minmod，`2` 为 van Leer。局部加密通过 `x_ref_lo/hi/split` 和 `r_ref_lo/hi/split` 配置，并以输出的实际网格面为准。
- `open_far_boundary=.true.` 把轴向远端设为开口。全域撞壁工况保留远端壁面，比较时区分贴壁气体压力与壁面动量通量。
