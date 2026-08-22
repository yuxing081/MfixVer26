# MFIX 26.1.2 JKR Validation

更新日期：2026-08-22

## 1. 验证对象

- MFIX 版本：26.1.2。
- JKR Git commit：`954ea0b9639592b9632c4f8fba92d23eabf42b99`。
- 验证框架父 commit：`5dd2fc9`；本阶段改动提交前运行证据为 dirty，完整 porcelain 状态保存在每次运行的 `metadata.json`。
- 串行 solver：`/home/yuxing/mfix/MfixVer26/code/build/mfixsolver`。
- 串行 solver SHA-256：`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`。
- 官方 stock serial module：`/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/mfixsolver.so`。
- 官方 stock module SHA-256：`3e95ef9efd32f8fe77ad991758adf16938bca7be617b7ec19852438d538be6b7`。
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
理论单测明确区分最大拉力点 `F=-F_c` 与位移控制临界分离端点 `F=-5F_c/9`；另外检查零重叠点 `F=-8F_c/9` 和零力平衡点。

### 2.2 输入生成与前置检查

所有共享物理和数值参数集中在 `validation/jkr_tools/config.py`。case 生成器在写入前检查：

- 必需关键字存在且只出现一次；
- JKR 工况的 `use_jkr` 与目的相符；stock Hertz 输入则完全移除 `use_jkr`、`jkr_*` 和表面能扩展关键字；
- 粒子数量、user scalar 列数正确；
- 粒子位于域内；
- 初始颗粒对重叠或墙面间隙符合配置；
- `dt <= T_H / 10000`；
- 双精度 particle output、速度、角速度和 particle ID 输出开关存在。

标准 case 提交到 `validation/cases`；运行器逐字节校验后复制到新的时间戳结果目录，因此不会在资产目录运行或覆盖历史证据。

普通 Hertz/JKR approach case 使用统一初始表面间距
`0.05R = 5e-8 m`：wall 粒心 y=1.05 μm，pair 粒心
y=8.975/11.025 μm。restart 工况为了预先建立 PFT/PFR/state 历史，
保留各自的初始接触或正重叠，不套用该间距。

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

状态含义：PASS 表示有本次 MFIX26 数值证据；DEFERRED 表示按当前范围明确推迟；PLANNED 表示尚待实现；BLOCKED 表示存在生产使用阻塞。

| 类别 | 工况 | 比较对象 | 主要指标 | 判据 | 状态 |
|---|---|---|---|---|---|
| Unit | JKR reference/case assets/plotting | 独立公式与生成器 | 关键点、临界速度、e=1 动态、34 个输入文件、四类图 | 67 个单元测试整体通过 | PASS |
| Unit | case generation | 配置与生成文本 | 唯一关键字、粒子数、几何、dt、输出 | 任一前置检查失败即 FAIL | PASS |
| Unit | I/O and comparison | 合成 CSV/monitor/VTP | 精度元数据、插值、global ID、pair ID | 重复/缺失/时间错位必须拒绝 | PASS |
| A 原生保护 | pair Hertz normal | 官方 26.1.2 vs custom JKR-off | x、v、omega、重叠、碰撞力、末态 | 共同 VTP 帧及双精度末态完全一致 | PASS |
| A 原生保护 | wall Hertz normal | 官方 26.1.2 vs custom JKR-off | 同上 | 共同 VTP 帧及双精度末态完全一致 | PASS |
| A 原生保护 | oblique Hertz friction | 官方 26.1.2 vs custom JKR-off | 切向运动、omega、碰撞力、末态 | 共同 VTP 帧及双精度末态完全一致 | PASS |
| A 原生保护 | wall/pair Hertz elastic matched speeds | 官方 26.1.2 vs custom JKR-off | e=1、无黏附、纯弹性回弹 | 轨迹/末态一致且 e_eff 误差不超过 5e-3 | PASS |
| B 墙面 | wall_normal_capture/rebound | e=1 JKR 理论与 MFIX26 | 0.95 捕获、1.05 断颈、总法向力、四类图 | 行为、关键点符号和逐点力误差通过 | PASS |
| B 墙面 | wall_normal_velocity_sweep | 固定速度扩展 | e_eff、重叠、时长 | 当前不搜索转变区间 | DEFERRED |
| B 墙面 | oblique stick/rebound | 理论趋势与 MFIX26 | PFT、切向力、omega、能量 | 库仑限制和历史连续 | PLANNED |
| B 墙面 | sliding/rolling | 独立历史与 MFIX26 | PFT/PFR、滚动力矩、临界位移 | 无意外清零且限制正确 | PLANNED |
| C 颗粒对 | pair_normal_capture/rebound | e=1 JKR 理论与 MFIX26 | 0.95 捕获、1.05 断颈、总法向力、四类图 | 行为、关键点符号和逐点力误差通过 | PASS |
| C 颗粒对 | pair_normal_velocity_sweep | 固定速度扩展 | e_eff、重叠、时长 | 当前不搜索转变区间 | DEFERRED |
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

`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T151407.503647Z-5944/unit-tests.log`

结果：67 tests，全部通过；覆盖 17 个 canonical case 的逐字节同步、15 个普通 approach case 的 0.05R 几何、e=1 理论捕获/回弹以及四类 contact PNG。

### 4.2 原生 Hertz 保护

命令：

    ./validation/run_validation.sh hertz

证据目录：

`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T135425.285145Z-2946`

对照双方使用逐字节相同的 stock-compatible `mfix.dat` 和
`particle_input.dat`。输入保留 `des_coll_model='HERTZIAN'`、材料、网格、
初始条件和输出设置，但不包含 stock MFIX26 不认识的 `use_jkr`、`jkr_*`
或表面能扩展关键字。

| 工况 | `e_eff` | 最大重叠 (m) | 接触时长 (s) | 最大轨迹误差 | 结果 |
|---|---:|---:|---:|---:|---|
| pair_hertz_normal | 0.59999436 | 1.7400680e-8 | 6.4998258e-8 | 0（x/v/omega/force） | PASS |
| wall_hertz_normal | 0.60566479 | 1.5191520e-8 | 5.5996958e-8 | 0（x/v/omega/force） | PASS |
| wall_hertz_oblique | 0.60566479 | 1.5191520e-8 | 5.5996958e-8 | 0（x/v/omega/force） | PASS |
| wall_hertz_elastic_capture | 1.00000468 | — | — | 0（x/v/omega/force） | PASS |
| wall_hertz_elastic_rebound | 1.00000061 | — | — | 0（x/v/omega/force） | PASS |
| pair_hertz_elastic_capture | 1.00000000 | — | — | 0（x/v/omega/force） | PASS |
| pair_hertz_elastic_rebound | 1.00000000 | — | — | 0（x/v/omega/force） | PASS |

迁移 solver 每个工况额外写出 `PARTICLES_0501.vtp`；其物理时间和内容均
与共同末帧 `PARTICLES_0500.vtp` 逐字节相同，因此记录为输出调度差异，
不当作物理差异。比较器只允许位于末尾且与最后共同帧完全相同的额外帧；
共同帧缺失、字段差异或非重复额外帧仍判为 FAIL。

官方 module SHA-256 为
`3e95ef9efd32f8fe77ad991758adf16938bca7be617b7ec19852438d538be6b7`；
迁移 solver SHA-256 为
`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`。

### 4.3 标准 case 资产与基础 JKR normal

维护命令：

    python3 -m validation.generate_cases

`validation/cases` 当前提交 17 个标准 case，每个目录仅含生成的
`mfix.dat` 与 `particle_input.dat`。67 个单元测试会将这 34 个文件与
当前生成器逐字节比较；运行器先执行同一同步检查，再复制到独立时间戳结果
目录，不在 case 资产目录内运行或覆盖历史证据。

15 个普通 approach case 的初始间距均由集中数值配置显式传入并经单元测试
逐一检查；wall/pair 均为 `0.05R`。两个 restart case 不属于 approach
case，保留用于建立历史的专用初始接触状态。


材料档案显式拆分为 `elastic_jkr`（e=1.0）与
`damped_frictional`（e=0.6）。e=1 JKR 使用固定 0.95/1.05 倍
adhesion-only critical speed；e=0.6 历史 normal case 已改名为
`*_normal_damped_capture/rebound`，保留 0.95/3.0 倍。本阶段按要求不做
自动临界搜索或自适应加密。

统一 `all` 证据目录：
`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T151407.503647Z-5944`，
unit、7 个 Hertz、8 个 JKR normal 和 2 个 restart case 共 18 项全部 PASS。

| e=1 工况 | 速度因子 | 预期/结果 | 力最大绝对误差 (N) | 最大误差时间 (s) |
|---|---:|---|---:|---:|
| wall_normal_capture | 0.95 | 捕获 / PASS | 1.078e-10 | 9.840e-7 |
| wall_normal_rebound | 1.05 | 断颈回弹 / PASS | 1.157e-10 | 1.640e-7 |
| pair_normal_capture | 0.95 | 捕获 / PASS | 5.445e-11 | 1.782e-6 |
| pair_normal_rebound | 1.05 | 断颈回弹 / PASS | 5.458e-11 | 1.980e-7 |

每个 JKR case 写出 `mfix_contact.csv`、
`python_reference_contact.csv` 和四张叠图：速度–时间、角速度–时间、
法向力–重叠、接触半径–时间。本轮 8 个 JKR case 共生成 16 份 CSV 和
32 张 PNG。接触半径的 MFIX 曲线由 MFIX 重叠按同一审查后的 JKR 几何关系
反算，图例明确标为 `MFIX26 overlap-inferred`；VTP 仍只作为轨迹与趋势
证据。metadata 记录材料、临界速度/因子、四个理论关键点、容差、首个差异、
最大误差时间以及 CSV/PNG 绝对路径。
四类图均保留首次接触前的 approach 时段；在法向力–重叠图中该段为零力，
在接触半径–时间图中该段明确绘为零半径，不再因 inactive/None 过滤掉。


### 4.4 串行 restart

最终证据目录：

`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T114550.104470Z-132230`

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

- 基础 normal 捕获/回弹、法向力、理论轨迹和绘图已执行；能量分解、异质表面能和 DMP 矩阵尚未完成，壁面几何后移。旧 MFIX23 或迁移期记录不会自动转为 PASS。
- 当前 restart 正向判据覆盖双精度末态，但尚未在 restart 边界导出并逐字段比较 NEIGHBORS、PFT、PFR、JKR_CONTACT_STATE、墙面 facet 和 wall history；这是进入完整 restart 验收前的主要缺口。
- Float32 VTP 不能用于边界历史的严格相等判断。
- 当前实现使用固定墙面历史槽；预期工况中单颗粒不会同时接触超过 8 个壁面。超过 8 个壁面接触暂不测试，记录为范围限制。
- 本阶段不执行自动临界搜索；基础 contact 的四类标准图已接入。能量分解和更宽固定速度矩阵仍待后续需求确认。
- DMP 首先只验收 2 ranks；更多进程数在 2-rank 通过后再考虑。
- SuperDEM/GSP 异形颗粒不在当前 JKR 范围内。

## 6. 总结

已完成并通过：

- 集中配置、17 个 Git 跟踪 canonical case、逐字节同步和运行时复制；
- 独立 JKR 理论参考、CSV/monitor/VTP 读取和 global-ID 比较基础；
- 不覆盖历史结果的统一运行入口与完整元数据；
- e=1 wall/pair 0.95× 捕获、1.05× 断颈回弹和四类 contact 图；
- 对应 e=1 Hertz 纯弹性回弹，以及保留并改名的 e=0.6 damped case；
- serial 颗粒对和斜碰墙面 RESTART_1 末态回归；
- RESTART_2 敏感性负对照；
- 官方 stock MFIX26 与迁移 solver 的 7 个 Hertz 保护对照；
- JKR 最大拉力、临界分离、零重叠和平衡关键点单测；
尚未通过且阻塞完整生产验收：

- 基础墙面/颗粒对能量分解（自动临界搜索按当前要求不执行）；
- 异质表面能；
- 2-rank DMP（随后进行）；
- restart 边界内部历史的双精度逐字段连续性。
- 壁面几何及多壁面同时接触（后移，不阻塞前两项）。

下一阶段先完成异质表面能，再进行 2-rank DMP 和 restart 双精度边界诊断；壁面几何后移。不以旧版本结果或单纯运行完成替代数值证据。
