# MFIX 26.1.2 JKR 下一阶段测试迁移 Prompt

继续完善：

`/home/yuxing/mfix/MfixVer26`

只读参考：

`/home/yuxing/mfix/P3/JKR/versions/test`

不得修改 `versions` 中的任何内容，也不得把 MFIX23 solver、results、
cache 或数值输出复制到 MFIX26 并作为真值。可以学习旧版的测试组织方式、
case 设计、比较流程和绘图逻辑。

当前基线 commit 为 `7b98d5f`。开始前先检查实际 Git 状态和最新 commit；
如果与本文不同，以只读检查结果为准并记录，不要重置或覆盖用户修改。

## 一、本阶段目标

按以下顺序推进，并在每个阶段通过后再进入下一阶段：

1. 重构新运行的结果目录，使其按测试类别和物理子类组织；
2. 建立并运行 serial 异质表面能测试；
3. 建立并运行 2-rank DMP 测试；
4. 设计并实现默认关闭的 restart 边界双精度诊断，完成 restart 内部历史比较；
5. 壁面多面同时接触和复杂 wall geometry 暂时后移，只预留目录与计划。

不要执行自动临界区间搜索。普通接近工况继续使用
`initial_gap = 0.05R = 5e-8 m`；restart/history 工况为建立有效历史，
允许使用明确记录的初始接触或正重叠。

## 二、结果目录契约

参考 MFIX23 的可浏览结构，结果直接保存在固定分类目录：

```text
validation/results/
├── RUN_INFO.txt
├── metadata.json
├── summary.csv
├── summary.json
├── unit/
│   └── unit-tests.log
├── contact_physics/
│   ├── pair_normal/<case_name>/<execution_mode>/
│   ├── pair_oblique/<case_name>/<execution_mode>/
│   ├── wall_normal/<case_name>/<execution_mode>/
│   └── wall_oblique/<case_name>/<execution_mode>/
├── dmp/
│   ├── pair_cross/<case_name>/{serial,dmp}/
│   ├── pair_new_contact/<case_name>/{serial,dmp}/
│   ├── active_contact_migration/<case_name>/{serial,dmp}/
│   ├── wall_cross/<case_name>/{serial,dmp}/
│   └── cloud/<case_name>/{serial,dmp}/
├── restart/
│   ├── pair_history/<case_name>/{continuous,stage1,restart1,restart2}/
│   └── wall_history/<case_name>/{continuous,stage1,restart1,restart2}/
└── wall_geometry/
    ├── single_face/<case_name>/serial/
    ├── multi_face/<case_name>/serial/
    └── facet_seam/<case_name>/serial/
```

目录层级语义固定为：

1. 第一层是测试类别，如 `contact_physics`、`dmp`、`restart`、
   `wall_geometry`；
2. 第二层是物理子类，如 `pair_normal`、`pair_oblique`、
   `wall_normal`、`wall_oblique`、`pair_cross`；
3. 第三层是明确的 case 名称；
4. 最内层是执行模式或比较对象，如 `serial`、`dmp`、
   `official`、`migrated`、`continuous`、`restart1`。

Hertz、JKR 和异质表面能由 case 名称及 metadata 明确区分，不在
`contact_physics` 下额外增加模型层，以便打开目录后直接看到正碰、斜碰等
物理子类。

不要为追求完全对称而创建空目录。只在实际运行时创建相应路径。每次调用
统一 runner 时先清理上一份结果快照，再写入一套自洽的新结果；Git commit、
canonical case 和 solver hash 用于重新生成需要的旧结果。
`RUN_INFO.txt` 记录运行时间、MFIX 版本、Git commit/dirty、solver hash、
suite、执行模式和 MPI 进程数。所有 summary 和 metadata 中的路径必须与
当前分类目录一致。

在集中配置或 case registry 中显式保存结果类别与子类，不要在多个 runner
分支中重复硬编码字符串。增加单元测试验证：

- 每个 canonical case 恰好映射到一个合法结果路径；
- 相同 case 与执行模式始终映射到相同相对路径；
- case 名和执行模式不得丢失；
- 两个不同 case 不得映射到同一目录；
- 清理上一快照时保留 `.gitkeep`，且新结果中不得混入旧运行残留；
- summary 中路径与实际目录一致。

## 三、serial 异质表面能

先只读审查 MFIX26 生产源码中以下契约，再确定 case：

- `DES_USR_VAR(1,:)` 的颗粒表面能语义；
- 颗粒值未设置、为零、为负或非法时的实际 fallback 行为；
- 颗粒–颗粒与颗粒–墙面的 `gamma_eff` 组合；
- global particle ID 是否只用于身份而不影响组合结果。

独立参考继续采用经源码确认后的显式函数。不得为了测试通过而改变公式。
至少建立：

1. `pair_gamma_homogeneous`：两个颗粒 gamma 相同，退化为现有均质结果；
2. `pair_gamma_low_high`：相同固相内 low/high gamma；
3. `pair_gamma_high_low`：交换输入行和 global ID 后的物理对称对照；
4. `pair_gamma_id_swap`：保持物理初态，仅交换 global ID；
5. `wall_gamma_heterogeneous`：颗粒 gamma 与墙面 gamma 不同；
6. `pair_gamma_fallback`：未设置值回退到 phase 默认值；
7. 仅在源码契约明确支持时增加 invalid-value case；若生产代码会显式拒绝，
   应把“正确报错”作为 PASS，而不是强行要求 fallback。

这些基础本构 case 使用 `elastic_jkr`、e=1 和 0.05R 初始间距。
每个 case 记录 `gamma_i`、`gamma_j/gamma_wall`、独立计算的
`gamma_eff`、`W=2 gamma_eff`、`a0`、`F_c`、关键位移和容差。

通过标准至少包括：

- 同 gamma 与均质 case 在规定容差内一致；
- low-high 与 high-low 的力–重叠、轨迹和断颈行为一致；
- 交换 global ID 后按物理粒子重新对齐，结果保持对称；
- MFIX 总法向力与独立参考曲线一致；
- fallback 或输入拒绝行为与生产源码契约完全一致；
- 输出首个差异、最大误差及发生时间，而不是只看 solver exit code。

## 四、2-rank DMP

serial 异质表面能全部通过后再开始。固定使用 2 ranks 和 2×1×1 分区，
不复制 MFIX23 的 dmp8 设置。

依次建立：

1. `pair_cross_jkr`：不同 rank 上的颗粒建立接触；
2. `pair_new_contact_jkr`：跨区后首次建立接触；
3. `pair_active_contact_migration`：已有非零物理 PFT/PFR 和 JKR state
   的接触迁移或发生 ownership reverse；
4. `wall_cross_jkr`：颗粒跨 rank 后继续保持墙面接触；
5. `cloud_collision_jkr`：小颗粒群 serial/DMP 集成对照。

不得用固定时间、固定 particle ID 或验证专用 history seed 修改生产物理。
需要非零 PFT/PFR 时，应通过斜碰、滚动或持续接触自然形成。

所有比较必须：

- 粒子按 global ID 对齐；
- pair contact 按 canonical `(min_gid,max_gid)` 对齐；
- wall contact 按 global particle ID 与稳定 facet/物理墙面 ID 对齐；
- 比较轨迹、速度、角速度、力、接触集合、PFT、PFR 和 JKR state；
- 检查接触不丢失、不重复、力不双算、ownership reverse 向量符号正确；
- 特别观察 `mpi_unpack_des_mod.f` 邻居 fallback 分支的 neighbor 定位；
- serial/DMP 超差时报告第一个物理差异，而不是局部数组索引差异。

## 五、restart 双精度边界诊断

完成 2-rank 基线后，为生产代码增加最小、默认关闭、无 case-specific
条件的诊断接口。不要复用 `DES_USR_VAR(2:)`，不要依赖 Float32 VTP。

诊断应以双精度文本或 stream 输出：

- time、phase tag、rank；
- owner/neigh global ID 或 wall facet ID；
- local slot 仅作辅助信息；
- overlap、normal/tangential force；
- PFT、PFR；
- JKR contact state；
- neighbor rebuild 前后、pending queue 恢复后、首次 force evaluation 后。

比较 continuous、split+RESTART_1 与 RESTART_2 敏感性。RESTART_1 必须在
相同物理时刻逐字段一致；RESTART_2 只作负对照，不要求保存 JKR history。
诊断开关默认 false，不得改变未启用时的数值轨迹。

## 六、壁面几何范围

`wall_geometry` 保留为一级类别，但本阶段不实现多壁面同时接触。
后续再做 single face、edge/corner、真实 facet seam，并检查接触重复、
合力方向、facet 切换及 wall PFT/PFR/state 连续性。单颗粒超过 8 个墙面
接触继续记录为当前范围外。

## 七、实施和验收

每完成一个阶段立即更新：

- `docs/TEST_PROGRESS.md`；
- `docs/VALIDATION.md`；
- `validation/results` 中的 `RUN_INFO.txt`、`metadata.json`、
  `summary.csv` 和 `summary.json`。

每阶段运行：

- `py_compile`；
- 全部单元测试；
- canonical case 字节同步；
- 相应最小真实工况；
- 阶段统一入口；
- 最后运行 `./validation/run_validation.sh all`。

失败的独立测试应继续执行并汇总；任一必需测试失败时统一入口返回非零。
大型原始输出不提交 Git，只提交 case、配置、生成器、比较器、小型 CSV/图、
摘要和文档。提交前运行 `git diff --check`，检查 Git 状态并审查所有
参数、日志、比较证据和图形。

