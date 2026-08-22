# MFIX 26.1.2 JKR Validation

更新日期：2026-08-22

## 1. 验证对象

- MFIX 版本：26.1.2。
- JKR Git commit：`954ea0b9639592b9632c4f8fba92d23eabf42b99`。
- 当前验证工作树：dirty；新增验证框架尚未提交，完整 porcelain 状态保存在每次运行的 `metadata.json`。
- 串行 solver：`/home/yuxing/mfix/MfixVer26/code/build/mfixsolver`。
- 串行 solver SHA-256：`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`。
- 构建编译器：`/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/gfortran`，GNU Fortran 15.3.0。
- 当前已验收执行环境：serial，1 process，分区 1×1×1。
- DMP 目标环境：2 ranks，分区 2×1×1；尚未在本体系中重新验收。

每次运行必须重新记录 MFIX 版本、Git commit、dirty 状态、solver 路径与 SHA-256、实际构建编译器、执行模式、MPI 进程数、分区和全部集中配置。大型结果只保留绝对路径和哈希，不纳入 Git。

## 2. 测试方法

### 2.1 参考层级

判断优先级固定为：

1. 解析理论或独立 Python 参考模型；
2. 未修改的 MFIX 26.1.2 原生 Hertz solver；
3. 自定义 solver 在 `use_jkr=.False.` 时的原生保护对照；
4. serial 与 DMP 对照；
5. continuous 与 RESTART_1 对照；
6. MFIX 23.4.1 JKR 结果仅作迁移回归参考，不作为理论真值。

当前独立参考明确采用以下 MFIX 生产契约：

- 异质表面能：`gamma_eff = sqrt(gamma_i * gamma_j)`；
- 颗粒–墙面：`gamma_eff = sqrt(gamma_particle * gamma_wall)`；
- 文献 work of adhesion：`W = 2 gamma_eff`；
- 颗粒对有效半径：`R_eff = R_i R_j / (R_i + R_j)`；
- 颗粒–墙面有效半径：`R_eff = R`；
- 平衡接触半径：`a0 = (9 pi gamma_eff R_eff^2 / E_eff)^(1/3)`；
- 临界分离量：`delta_c = a0^2 / (2 6^(1/3) R_eff)`；
- 拉脱力幅值：`F_c = 3 pi R_eff gamma_eff`；
- 正向为压缩，负向为拉伸。

这些公式从旧 `versions/test/jkr_tools/reference.py` 审查后重写为显式标量函数；没有直接复制旧运行结果。

### 2.2 输入生成与前置检查

所有共享物理和数值参数集中在 `validation/jkr_tools/config.py`。case 生成器在写入前检查：

- 必需关键字存在且只出现一次；
- `use_jkr` 与工况目的相符；
- 粒子数量、user scalar 列数正确；
- 粒子位于域内；
- 初始颗粒对重叠或墙面间隙符合配置；
- `dt <= T_H / 10000`；
- 双精度 particle output、速度、角速度和 particle ID 输出开关存在。

写入目标目录必须不存在，因此生成器不会覆盖历史证据。

### 2.3 比较与容差原则

比较器返回 PASS/FAIL、首个超差字段/索引/物理时间，以及每个字段的最大绝对误差和发生位置。DMP 粒子按 global particle ID 对齐；接触按 `(min_gid, max_gid)` 物理 pair ID 对齐，重复、缺失或时间不一致直接失败。

初始容差分层如下；在真实基准建立前不得为通过测试而放宽：

| 数据源 | 位置绝对/相对 | 速度绝对/相对 | 角速度绝对/相对 |
|---|---:|---:|---:|
| 双精度内部/CSV | 1e-13 / 1e-10 | 1e-10 / 1e-9 | 1e-4 / 1e-9 |
| monitor 文本 | 1e-11 / 1e-7 | 1e-7 / 1e-6 | 按字段单独设定 |
| Float32 VTP | 2e-11 / 1e-7 | 1e-4 / 1e-5 | 4 / 1e-5 |

Float32 VTP 只用于轨迹和可视化，不能替代 restart 边界的双精度诊断。

## 3. 测试矩阵

状态含义：PASS 表示有本次 MFIX26 数值证据；READY 表示工具已实现但尚未运行完整工况；PLANNED 表示尚待实现；BLOCKED 表示存在生产使用阻塞。

| 类别 | 工况 | 比较对象 | 主要指标 | 判据 | 状态 |
|---|---|---|---|---|---|
| Unit | JKR reference | 独立公式与审查基准 | gamma、W、R/E/G、a0、delta_c、Fc、临界速度 | 40 个单元测试整体通过 | PASS |
| Unit | case generation | 配置与生成文本 | 唯一关键字、粒子数、几何、dt、输出 | 任一前置检查失败即 FAIL | PASS |
| Unit | I/O and comparison | 合成 CSV/monitor/VTP | 精度元数据、插值、global ID、pair ID | 重复/缺失/时间错位必须拒绝 | PASS |
| A 原生保护 | pair Hertz normal | 官方 26.1.2 vs custom JKR-off | x、v、omega、接触时间、重叠、力 | 双精度分层容差内一致 | PLANNED |
| A 原生保护 | wall Hertz normal | 官方 26.1.2 vs custom JKR-off | 同上 | 双精度分层容差内一致 | PLANNED |
| A 原生保护 | oblique Hertz friction | 官方 26.1.2 vs custom JKR-off | 切向冲量、omega、切向力 | 双精度分层容差内一致 | PLANNED |
| B 墙面 | wall_normal_capture/rebound | JKR 理论与 MFIX26 | 捕获、回弹、重叠、时长、Fc、能量 | 行为正确且误差不超容差 | PLANNED |
| B 墙面 | wall_normal_velocity_sweep | 自适应速度扫描 | 临界区间、e_eff、重叠、时长 | 覆盖捕获/临界/回弹并加密一次 | PLANNED |
| B 墙面 | oblique stick/rebound | 理论趋势与 MFIX26 | PFT、切向力、omega、能量 | 库仑限制和历史连续 | PLANNED |
| B 墙面 | sliding/rolling | 独立历史与 MFIX26 | PFT/PFR、滚动力矩、临界位移 | 无意外清零且限制正确 | PLANNED |
| C 颗粒对 | pair_normal_capture/rebound | JKR 理论与 MFIX26 | 捕获、回弹、重叠、时长、Fc | 行为正确且误差不超容差 | PLANNED |
| C 颗粒对 | pair_normal_velocity_sweep | 自适应速度扫描 | 临界区间、e_eff、重叠、时长 | 覆盖三类行为并加密一次 | PLANNED |
| C 颗粒对 | oblique/sliding/rolling | 独立历史与 MFIX26 | PFT/PFR、力、力矩、omega | 历史和符号连续 | PLANNED |
| C 颗粒对 | tensile_neck_break | JKR 状态机与理论 delta_c | 压缩、卸载、负重叠、断颈 | 完整经历六阶段且断裂位置正确 | PLANNED |
| D 异质表面能 | homogeneous/heterogeneous/swap/fallback | 独立 gamma_eff 函数 | 力曲线、交换对称、默认回退 | global ID 交换后物理对称 | PLANNED |
| E 壁面几何 | single_face/edge/corner/seam | 独立壁面贡献、serial/DMP | 接触数、方向、facet、PFT/PFR | 不重复、不丢失、合力正确 | PLANNED |
| F DMP | pair cross/new/migration | serial vs 2-rank | global ID 轨迹、pair state、PFT/PFR | 不丢失/重复/双算且容差内一致 | PLANNED |
| F DMP | wall cross | serial vs 2-rank | facet、wall state、PFT/PFR | rank 迁移后连续 | PLANNED |
| F DMP | cloud collision | serial vs 2-rank | 轨迹、能量、接触集合 | global ID 对齐后容差内一致 | PLANNED |
| G Restart | restart_pair_pfr | continuous vs RESTART_1 | 双精度 particle output | 逐字节一致 | PASS |
| G Restart | restart_wall_oblique | continuous vs RESTART_1 | 双精度 particle output | 逐字节一致 | PASS |
| G Restart | RESTART_2 sensitivity | continuous vs RESTART_2 | 双精度 particle output | 差异被明确记录，不作正向通过条件 | PASS |
| G Restart | internal boundary history | 连续/重启边界双精度诊断 | pair ID、neighbors、PFT/PFR/state、facet | 首次重建/力计算前后逐字段一致 | PLANNED |

## 4. 结果

### 4.1 单元测试

命令：

    ./validation/run_validation.sh unit

结果目录：

`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T095448.025415Z-107484`

结果：40 tests，全部通过。覆盖 8 个 case-generation、9 个 comparison、5 个 I/O、15 个 reference 和 3 个 plotting 数据检查。

### 4.2 串行 restart

最终证据目录：

`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T095448.025415Z-107484`

运行信息：

- MFIX 26.1.2；
- solver SHA-256：`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`；
- GNU Fortran 15.3.0，来源由 solver 构建目录的 CMakeCache 解析；
- serial，1 process，1×1×1；
- continuous、stage1、RESTART_1、RESTART_2 共 8 次 solver 执行，退出码均为 0。

关键结果：

| 工况 | continuous SHA-256 | RESTART_1 | RESTART_2 敏感性 |
|---|---|---|---|
| restart_pair_pfr | `b80cf63a…454806d` | 相同，PASS | `3142dea0…46e980`，不同 |
| restart_wall_oblique | `b3ad2958…0f402` | 相同，PASS | `2c9e1cb9…941d8d2`，不同 |

本阶段的数值判据不是“solver 能运行”，而是 continuous 与 RESTART_1 的双精度末态文件逐字节一致。RESTART_2 的不同说明负对照确实对历史数据敏感。

## 5. 已知限制

- 尚未重新执行原生 Hertz 保护、基础捕获/回弹、速度扫描、异质表面能、壁面几何和 DMP 矩阵；旧 MFIX23 或迁移期记录不会自动转为 PASS。
- 当前 restart 正向判据覆盖双精度末态，但尚未在 restart 边界导出并逐字段比较 NEIGHBORS、PFT、PFR、JKR_CONTACT_STATE、墙面 facet 和 wall history；这是进入完整 restart 验收前的主要缺口。
- Float32 VTP 不能用于边界历史的严格相等判断。
- 当前实现使用固定墙面历史槽；预期工况中单颗粒不会同时接触超过 8 个壁面。超过 8 个壁面接触暂不测试，记录为范围限制。
- 速度扫描 CSV、自动临界加密和五类标准图尚未接入统一入口。
- DMP 首先只验收 2 ranks；更多进程数在 2-rank 通过后再考虑。
- SuperDEM/GSP 异形颗粒不在当前 JKR 范围内。

## 6. 总结

已完成并通过：

- 集中配置、运行时 case 生成和输入前置检查；
- 独立 JKR 理论参考、CSV/monitor/VTP 读取和 global-ID 比较基础；
- 不覆盖历史结果的统一运行入口与完整元数据；
- serial 颗粒对和斜碰墙面 RESTART_1 末态回归；
- RESTART_2 敏感性负对照。

尚未通过且阻塞完整生产验收：

- MFIX 26.1.2 官方 Hertz 保护对照；
- 基础墙面/颗粒对物理与速度扫描；
- 异质表面能和壁面几何；
- 2-rank DMP；
- restart 边界内部历史的双精度逐字段连续性。

下一阶段按上述顺序实现和运行，不以旧版本结果或单纯运行完成替代数值证据。
