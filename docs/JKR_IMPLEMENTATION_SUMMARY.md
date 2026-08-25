# MFIX 26.1.2 球形颗粒 Chen/JKR 实现汇总

更新日期：2026-08-25
代码基线：`fdf4c2af3f0085f5b5428a7b34a6cdb2304dc9b0`
验证基线：`aedb8f9c0487d4a7edbf129d45981c126921ee85`

## 1. 结论与范围

本项目以 MFIX 26.1.2 官方球形 Hertzian DEM 为主体，将 MFIX 23.4.1
分支中的 JKR 接触物理迁移并按新版数据结构重写。当前实现覆盖：

- 球形颗粒–颗粒和颗粒–墙面的 JKR 法向接触、负重叠 neck、黏附滑动上限；
- Chen/Marshall 滚动模型的弹性 spring-slider 部分：滚动位移、临界滚动角
  和滚动力矩；
- 逐相、逐颗粒及颗粒–墙面异质表面能；
- MFIX 26.1.2 邻居重建、粒子排序、STL facet、DMP owner migration 和
  `RESTART_1` 所需的接触身份与历史连续性。

当前证据支持“在上述范围内较完整地复现了 Chen 2019 采用的球形颗粒
JKR 接触和滚动模型的弹性 spring-slider 部分”，不等于实现了 Chen 2019
的全部 fast-DEM 方法。
具体而言，生产代码不自动执行论文 Eq. (10) 的降刚度表面能缩放，没有为
Eq. (11)/(12) 单独保存原始表面能 `gamma_O`，也没有实现 twisting torque
Eq. (1c)/(3b) 或复现论文中的宏观堆积统计。

JKR 仅支持普通球形 Hertzian DEM。`SuperDEM`、显式/隐式 GSP、LSD 和
原生 `USE_COHESION` 不进入 JKR 分支；`USE_JKR=.FALSE.` 是保持 MFIX
26.1.2 原生行为的兼容边界。

## 2. 基线与迁移原则

| 对照对象 | 本项目的处理 |
|---|---|
| MFIX 23.4.1 JKR，commit `ea55027498a4b939a2982cbb77d137201b034ffc` | 保留经审计的球形 JKR、PFT/PFR 和滚动物理；不复制旧版二维邻居身份、`NEIGHBOR_INDEX` restart 和已知缺陷 |
| MFIX 26.1.2 官方源码 | 作为生产结构基线；复用其 global particle ID、contact slot、STL 过滤、owner/ghost、pack/unpack 和 restart 前缀 |
| 本仓库基线 commit `2cb0da5` | 相对该 commit 的 `code/` 差异定义本次迁移的实际源码增量 |

接触身份遵守两个约束：一个物理接触只贡献一次力；同一物理接触只有一份
权威历史。颗粒对以 canonical global-ID pair 识别，墙接触复用 MFIX
`WALL_COLLISION_FACET_ID` slot，并只在已证明原生 facet identity 无法保持
连续性的共面 seam 上最小扩展。

`code/read_namelist.f` 与 `model/read_namelist.f` 逐字节相同，只是自定义构建
overlay 所需的官方文件，不属于 JKR 语义修改。

## 3. 输入、派生量和历史变量

### 3.1 新增输入

声明和默认值位于 `code/desnamelist.inc:50-53` 与
`code/des_init_namelist.f:1494-1533`，合法组合和值域检查位于
`code/check_solids_dem.f:500-516,921-963`。

| 输入 | 默认值 | 含义 |
|---|---:|---|
| `USE_JKR` | `.FALSE.` | 开启球形颗粒 JKR 接触和 Chen 滚动 |
| `DES_SURFACE_ENERGY(M)` | `UNDEFINED` | 固相 `M` 的表面能 |
| `DES_SURFACE_ENERGY_WALL` | `UNDEFINED` | 墙面表面能 |
| `JKR_THETAR_CRIT` | `1.0e-2` | 颗粒对临界滚动角 |
| `WALL_JKR_THETAR_CRIT` | `1.0e-2` | 颗粒–墙面临界滚动角 |
| `JKR_ALPHA_SAFETY` | `1.0` | 拉伸 cutoff 相对 `delta_c` 的安全系数，允许范围 0.7–1.0 |

`DES_USR_VAR(1,particle)` 是可选的逐颗粒表面能；正值覆盖所属相的默认值，
非正值回退并规范化为相默认值。访问和回退分别见
`code/discretelement_mod.f:1090-1103` 与 `code/calc_force_dem.f:130-138`。

### 3.2 派生系数

`code/discretelement_mod.f:380-398` 新增以下派生量：

- 表面能组合、平衡接触半径、拉脱力和临界分离量：`JKR_GAMAN/AoN/FcN/DeltacN`
  及墙面对应数组；
- 以接触半径为运行时尺度的 `HERT_KN_JKR/HERT_KT_JKR` 及墙面数组；
- 按 `sqrt(a)` 缩放的 `DES_ETAN_JKR/DES_ETAT_JKR` 及墙面数组。

相级默认系数在 `code/check_solids_dem.f:1130-1155,1241-1255` 初始化。
实际接触允许逐颗粒表面能，因此 `a0`、`delta_c` 和 `Fc` 会在力循环中按
实际粒子重新计算，而不是无条件使用相级矩阵。

### 3.3 接触历史

| 历史 | 数据布局 | 方向契约 |
|---|---|---|
| `PFT_NEIGHBOR(:,CC)` | MFIX 原生颗粒对切向位移 | owner-directed；pair owner 反转时取反 |
| `PFR_NEIGHBOR(:,CC)` | 新增颗粒对滚动位移 | pair-invariant；owner 反转时不变 |
| `JKR_CONTACT_STATE(CC)` | 新增颗粒对 neck 状态 | 无方向；owner 反转时不变 |
| `WALL_COLLISION_PFR(:,slot,particle)` | 新增墙面滚动位移 | 依附 MFIX facet slot |
| `JKR_WALL_CONTACT_STATE(slot,particle)` | 新增墙面 neck 状态 | 依附 MFIX facet slot |

声明位于 `code/discretelement_mod.f:254-259,712-714`。这些数组复用 MFIX
已有的 `CC`、particle row 和 wall slot 生命周期，没有另建平行的接触表。

### 3.4 新增或专门扩展的函数

| 函数 | 位置 | 职责 |
|---|---|---|
| `JKR_DAMPING_ALPHA` | `code/check_solids_dem.f:27-40` | 由 restitution 计算 Marshall JKR 阻尼系数 |
| `GET_DES_SURFACE_ENERGY` | `code/discretelement_mod.f:1090-1103` | 返回逐颗粒 override 或相默认表面能 |
| `NEWTON_RAPHSON_METHOD` | `code/calc_force_dem.f:781-814`；`code/calc_collision_wall_mod.f:897-930` | 求 Chen Eq. (6) 的稳定接触半径分支 |
| `EXCHANGE_NEIGHBOR_DATA_JKR` | `code/des_functions.f:622-680` | 球形 JKR 粒子排序时同步搬运双精度 PFT/PFR/state |
| `QUEUE_PAIR_HISTORY`、`GROW_PENDING_PAIR_HISTORY`、`RESTORE_PENDING_PAIR_HISTORY` | `code/mpi_unpack_des_mod.f:834-947` | 按 global-ID pair 暂存并在邻居重建后恢复 DMP/restart 历史 |
| `FACETS_SHARE_IGNORED_EDGE`、`RETAG_JKR_HISTORY_ACROSS_IGNORED_EDGE`、`RETAG_JKR_HISTORY_TO_CURRENT_FACET` | `code/calc_collision_wall_mod.f:932-1106` | 验证共面共享边并在唯一安全目标上重标记现有 MFIX wall slot |
| `JKR_DIAGNOSTICS_ACTIVE`、`JKR_DIAG_PAIR_EVENT` | `code/mpi_unpack_des_mod.f:949-1000` | 默认关闭的双精度验证事件流 |

## 4. Chen/JKR 公式与生产代码

以下采用代码约定：压缩 overlap `delta` 为正，`n` 从 owner 指向 neighbor，
`x=a/a0`。异质表面能采用

```text
gamma_eff = sqrt(gamma_i gamma_j),
gamma_eff = sqrt(gamma_particle gamma_wall)  (particle-wall).
```

独立参考将论文的 work of adhesion 写成 `W=2 gamma_eff`。颗粒对
`R_eff=Ri Rj/(Ri+Rj)`，颗粒–墙面 `R_eff=R`。

### 4.1 接触参数和接触半径

根据 Chen 2019 §2.1 和 Eq. (6)：

```text
a0      = (9 pi gamma_eff R_eff^2 / E_eff)^(1/3)
Fc      = 3 pi R_eff gamma_eff
delta_c = a0^2 / (2 6^(1/3) R_eff)
delta / (6^(1/3) delta_c) = 2 x^2 - (4/3) sqrt(x)
```

颗粒对参数见 `code/calc_force_dem.f:344-360`，墙面参数见
`code/calc_collision_wall_mod.f:256-277`。稳定接触半径分支由 Newton 迭代
求解，见 `code/calc_force_dem.f:781-814` 和
`code/calc_collision_wall_mod.f:897-930`。

正 overlap 建立 neck；负 overlap 只延续已有 state。分离超过
`JKR_ALPHA_SAFETY*delta_c` 时同时清除 PFT、PFR 和 state，见
`code/calc_force_dem.f:362-427`。墙面遵守相同规则，见
`code/calc_collision_wall_mod.f:602-618,674-683`。

### 4.2 法向力、刚度和阻尼

Chen Eq. (1a) 的生产形式为

```text
F_N = -4 Fc (x^3 - x^(3/2)) n - eta_N (v_rel . n) n.
```

法向和切向刚度对应 Chen §2.1：

```text
k_N = (4/3) E_eff a
k_T = 8 G_eff a.
```

生产实现见 `code/calc_force_dem.f:434-468`，墙面实现见
`code/calc_collision_wall_mod.f:609-683`。JKR 阻尼基系数乘 `sqrt(a)`；
restitution 到无量纲阻尼系数的多项式来自 Marshall 2009 Eq. (53)，见
`code/check_solids_dem.f:27-40`，不是 Chen 2019 新提出的公式。

### 4.3 切向力和滑动上限

PFT 继续使用 MFIX 26.1.2 的切平面位移更新。弹簧–阻尼切向力为

```text
F_T = -k_T PFT - eta_T v_T.
```

根据 Chen Eq. (1b)/(3a)，JKR Coulomb 上限写为

```text
|F_T| <= mu |4 Fc (x^3 - x^(3/2)) + 2 Fc|.
```

滑动时同时把 PFT 回写为与截断后切向力一致的值。颗粒对实现见
`code/calc_force_dem.f:470-509`，墙面实现见
`code/calc_collision_wall_mod.f:689-725`。

### 4.4 Chen 滚动

滚动历史 `PFR` 在接触切平面中累计并限制为

```text
|PFR| <= theta_crit R_eff.
```

滚动刚度和力矩为

```text
k_R = 4 Fc x^(3/2)
M_R = k_R (n x PFR),
```

因此饱和时得到 Chen Eq. (3c)：

```text
|M_R,crit| = 4 Fc x^(3/2) theta_crit R_eff.
```

颗粒对实现见 `code/calc_force_dem.f:652-684`，墙面实现见
`code/calc_collision_wall_mod.f:795-815`。颗粒对 `PFR` 与 global ID
编号无关，作用–反作用力矩相反。当前实现对应 Eq. (1d) 的弹性
spring-slider 部分，没有单独的显式滚动阻尼 `eta_R` 项。

## 5. MFIX 26.1.2 接触生命周期扩展

本节内容用于保证上述公式在邻居变化、并行迁移和 restart 后仍作用于同一
物理接触，本身不是 Chen 公式。

| 生命周期 | MFIX 26.1.2 机制及最小 JKR 扩展 | 主要位置 |
|---|---|---|
| allocate/init/grow/deallocate | PFR/state 跟随原 PFT、particle row 和 wall slot | `code/des_allocate_mod.f:270-294,607-719,737-908`；`code/des_init_arrays.f:128-135`；`code/allocate_arrays.f:2049-2064` |
| particle swap/sort | 同步交换 wall history；JKR 球形排序使用双精度 PFT/PFR/state 搬运，JKR-off 保留官方例程 | `code/des_functions.f:177-337,622-680,778-928` |
| neighbor rebuild | 以 global-ID pair 快照和恢复；反向恢复仅 PFT 变号 | `code/neighbour.f:21-202` |
| wall facet slot | 在 MFIX `GET_COLLISION/REMOVE_COLLISION` slot 内同步清零或保留 PFT/PFR/state | `code/calc_collision_wall_mod.f:602-604,1113-1280` |
| coplanar STL seam | 活动 JKR contact 去重；跨 ignored shared edge 时重标记原 slot 的 facet ID，不复制历史 | `code/calc_collision_wall_mod.f:465-537,932-1106` |
| DMP owner migration | wall history 跟随 particle packet；pair 以 global IDs、cells、PFT/PFR/state 进入 pending queue，重建后恢复 | `code/mpi_init_des_mod.f:162-170`；`code/mpi_pack_des_mod.f:415-421,527-531`；`code/mpi_unpack_des_mod.f:531-537,618-635,834-947` |
| restart | JKR 使用 DES restart 2.0，在完整官方 1.2 前缀后追加 pair global IDs、PFR/state 和 wall facet/PFT/PFR/state | `code/write_res0_des.f:253-283`；`code/read_res0_des.f:67-72,292-348` |
| neighbor coarse-grid 判据 | `USE_JKR` 时 half-cell 判据系数由官方 2 修正为 4；关闭 JKR 时保持 2 | `code/desgrid_mod.f:1116-1126` |

默认关闭的双精度诊断由环境变量 `MFIX_JKR_DIAGNOSTICS` 启用，定义在
`code/mpi_unpack_des_mod.f:949-1000`。它记录 global pair/facet identity、
力、PFT、PFR 和 state，只用于验证，不参与物理计算。

## 6. 验证 case 与证据

统一入口为 `validation/run_validation.sh`。标准输入位于 `validation/cases/`，
独立理论和比较实现位于 `validation/jkr_tools/`，完整方法、容差和逐 case
结果见 [`validation/docs/VALIDATION.md`](../validation/docs/VALIDATION.md)。

### 6.1 证据链

| 证据层 | case/位置 | 已证明的内容 | 不能外推的内容 |
|---|---|---|---|
| 独立公式单测 | `validation/tests/test_reference.py` | `a0`、`Fc`、`delta_c`、关键力点、稳定分支、gamma 组合和阻尼拟合 | MFIX 接触生命周期 |
| 原生保护 | `validation/cases/native_hertz/` | 官方 MFIX 26.1.2 与自定义 solver 的 JKR-off 轨迹和双精度末态一致 | JKR 公式本身正确 |
| 法向 JKR | `particle_pair/pair_normal_*`、`particle_wall/wall_normal_*` | 捕获/回弹、压缩、负 overlap neck、临界分离和阻尼轨迹 | 多接触宏观结构 |
| 异质表面能 | `pair_gamma_*`、`wall_gamma_heterogeneous` | 几何平均、交换对称、逐颗粒 override 和 fallback | 任意材料混合律 |
| 斜碰和滚动 | `pair_oblique_sliding_rolling`、`wall_oblique_sliding_rolling` | stick/slip、PFT/PFR、Coulomb 上限、滚动饱和、力矩和角速度响应 | twisting torque |
| true-STL | `wall_geometry/`、`wall_seam_history/` | face/edge/corner 贡献数、共面去重、双向跨 seam PFT/PFR/state 连续 | 任意复杂 CAD 拓扑 |
| DMP | `dmp/`、`cloud_8rank/`、`rolling_id_invariance/` | global-ID 接触身份、owner migration、2-rank 和指定 8-rank 串并行一致性 | 任意 rank 数和分解 |
| restart | `restart/` 和 `dmp_restart` suite | 串行 RESTART_1 逐字节连续；2-rank 边界 PFT/PFR/state 精确恢复并满足固定轨迹判据 | RESTART_2 历史连续性 |
| 集成 smoke | `aggregate_1000/initial_aggregate_1000` | 1000 粒子固定 8-rank 路径可运行，修复后串并行角速度和碰撞力一致 | Chen 论文 packing fraction/配位数复现 |

JKR 独立标量参考位于 `validation/jkr_tools/reference.py`，切向和滚动参考位于
`validation/jkr_tools/tangential.py`；它们不调用生产 Fortran 公式，避免用
同一实现验证自身。DMP/restart 使用 global particle ID 或 canonical pair ID
对齐，Float32 VTP 只承担轨迹证据，历史边界由双精度诊断比较。

### 6.2 当前干净回归

2026-08-25 回归执行时，production 和 validation 工作树均为 clean，基线为
production commit `fdf4c2af3f0085f5b5428a7b34a6cdb2304dc9b0` 和 validation
commit `aedb8f9c0487d4a7edbf129d45981c126921ee85`。执行命令：

```text
cd validation
./run_validation.sh all
```

结果为 47/47 case PASS。unit 共运行 141 项，其中 137 PASS，4 项因缺少
可选 `fparser` 依赖跳过；跳过项是 impact graph 生成/解析检查，不参与
JKR 数值验收。执行覆盖 serial、2 ranks/2x1x1 和 8 ranks/2x2x2。

证据快照：

- `validation/results/RUN_INFO.txt`：两仓库 commit 和 `dirty=False`；
- `validation/results/metadata.json`：serial solver SHA-256
  `9bb2fdc48ae178d336fdc19af463493b6322571f886bdc0029bad67cdf466f07`，
  DMP solver SHA-256
  `fc43ab674f99e150d955905de201f3aa0dd743bd19efe62a62c0d5fb5d5c9912`；
- `validation/results/summary.csv`：47 项均为 `True`；
- `validation/results/unit/unit-tests.log`：141 项 unit 的逐项记录。

## 7. 结论边界和后续扩展

当前实现与验证可支持球形颗粒 Chen/JKR 在既定 serial、DMP、STL 和
RESTART_1 范围内的生产使用，但以下内容仍需单独设计和验证：

- 自动应用 Chen Eq. (10) 的 `gamma_R=(E_R/E_O)^(2/5) gamma_O`；
- 在降刚度条件下按 Chen Eq. (11)/(12) 使用原始 `gamma_O` 计算滚动和滑动
  极限；当前模型只有一套生产表面能，不应把这两式写成已经实现；
- twisting torque、SuperDEM/GSP、超过默认 8 个同时 wall facet slots；
- 不同 MPI 分解和 owner 方向的系统扫描、长期能量审计；
- Chen 2019 论文的 packing fraction、coordination number 和统计分布复现。

因此，准确表述是：本项目完成了 Chen 2019 所采用的核心球形 JKR 法向、
黏附滑动和滚动弹性 spring-slider 接触关系，并在 MFIX 26.1.2 中建立了相应
历史、几何、DMP 和 restart 的连续性证据；尚未声称完整复现 Chen 2019 的
fast-DEM 参数缩放及宏观堆积研究。
