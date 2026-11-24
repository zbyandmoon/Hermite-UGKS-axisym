# 一维求解器

## 编译与运行

在 Linux 或 WSL 的 Bash 中，从 `1d` 目录运行。

```bash
make
./ugks1d cases/sod_kn001_gh34.namelist
```

执行 `./ugks1d` 时，程序读取当前目录的 `input.namelist`。也可传入 `cases` 下的其他 namelist 或自定义文件。

## 输入参数

| 参数 | 用法 |
|---|---|
| `case_id` | 2 为静水平衡，3 为 Sod 激波管，4 为 Fourier 热传导 |
| `nx`、`xlength` | 空间网格数与域长 |
| `nv`、`vmax`、`Tref` | 速度节点数、速度范围与参考温度 |
| `quad_type` | `newton_cotes` 或 `gauss_hermite` |
| `force_type` | `finite_difference` 或 `hermite`。后者用 `Nm` 设置展开阶数 |
| `knudsen`、`phi` | 参考 Knudsen 数与加速度 |
| `cfl`、`maxtime` | CFL 数与终止时间 |
| `TwL`、`TwR` | Fourier 算例的左右壁温 |

程序在当前运行目录写入 `macro*.dat`。字段和对应的输出时刻见文件表头及运行日志。
