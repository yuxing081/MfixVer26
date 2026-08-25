# MFIX 原生实现的已知不足与本项目处理

更新：2026-08-24

本文记录迁移和验证过程中由源码审查、定量工况或 MFIX 官方支持答复确认的
MFIX 原生行为。每节单独标明证据等级，不把本项目新增模型的实现错误归因于
MFIX 上游。这里的“原生”指本仓库 `model` 中的 MFIX-26.1.2 基线；`code` 是本项目的
迁移实现。行号对应本文更新时的源码。项目坚持优先复用 MFIX 原有数据结构，
只在缺陷已被独立工况复现且修复边界明确时增加行为。

## 1. 共面 STL facet 切换会重置墙面切向历史

### 1.1 错误发生的物理过程

一个连续平面常由多个共面 STL 三角形组成。颗粒沿墙面滑动时，物理接触
并未中断，但最近三角面会从 facet A 切换为 facet B。MFIX 用整数
`FACET_ID` 作为墙面接触历史的身份：

1. A 接触期间，切向相对位移被积分到 A 对应槽的 `WALL_COLLISION_PFT`。
2. 颗粒越过内部共面边后，B 的 `FACET_ID` 无法匹配 A 的槽。
3. MFIX 为 B 选择空闲槽，并把 PFT 初始化为零。因此静摩擦弹簧在物理接触
   未结束时被重新开始，切向力和力矩可能产生不连续。

`STL_NB_ANGLE` 相关逻辑只忽略低曲率相邻 facet 的内部边碰撞，并不迁移
历史。原生代码在 `model/des/calc_collision_wall_mod.f:423-427` 跳过内部边；
`GET_COLLISION` 在同文件 `:717-738` 仅按 `FACET_ID` 查找，在 `:763-765`
为新 facet 写入 ID 并把 PFT 清零。`REMOVE_COLLISION` 也只按 facet ID
操作，见 `:843-854`。

原生 MFIX 因而会丢失 PFT；23.4.1 上的旧 JKR 修改继续沿用该槽身份，并把
PFR 和 JKR neck state 也放入同一 facet 槽，所以问题扩大为 PFT、PFR 和
tensile-neck state 同时不连续。旧文件
`/home/yuxing/mfix/P3/JKR/versions/code/calc_collision_wall_mod.f:860-898`
按 facet ID 查找并清零三个状态，`:997-1010` 删除旧槽时也不迁移历史。

### 1.2 本项目如何处理

JKR 路径继续使用 MFIX 的 `WALL_COLLISION_FACET_ID` 唯一权威槽，没有新增
另一套 history identity：

- `code/calc_collision_wall_mod.f:465-507` 在 MFIX 删除旧 ignored-edge
  facet 前处理两个可能的 facet 遍历顺序。
- 同文件 `:938-992` 确认真正共享的 ignored edge；`:1000-1060` 和
  `:1068-1106` 要求唯一目标、相同 BC/STL group、同向法线及活动 JKR
  state，然后只重写原槽的 facet ID。PFT、PFR 和 neck state 留在原槽。
- 精确落在共面 shared edge 上的重复 JKR 力由同文件 `:509-537` 去重。

修复目前有意限制在球形 JKR 接触；原生 Hertz/LSD 路径保持 MFIX 行为。
新增 `validation/cases/wall_seam_history/wall_stl_seam_hertz_pft_reset`
负对照使用 `USE_JKR=.FALSE.`、Hertz 接触和非零切向速度。实测 facet 1→2
前的 `|PFT|=3.1030462688305876e-9 m`，切换后的 PFT 为严格零，保留比例为
0。定量判据位于 `validation/jkr_tools/seam_history.py:134-182`，工况定义位于
`validation/jkr_tools/config.py:434-445`。该 case 的 PASS 表示稳定复现已知
原生缺陷，不表示 Hertz seam 历史在物理上连续。

## 2. 半网格邻居跳过阈值把两倍直径实现成两倍半径

### 2.1 错误发生的物理过程

邻居搜索会根据颗粒位于 DES 网格单元的哪一半，跳过远侧的一层相邻单元。
设单元宽度为 `W`、最大颗粒半径为 `R`、邻居搜索系数为 `F`。远侧单元
能够安全跳过的条件是半个单元宽度超过搜索范围：

    W/2 > 2*R*F

即 `W > 4*R*F = 2*diameter*F`。原生注释也写明 DES 网格不应小于
`2*dia*rlm_factor`，见 `model/des/desgrid_mod.f:1115`；但紧随其后的实现
使用 `2*MAX_RADIUS*FACTOR_RLM`，见 `:1116-1118`。这实际在
`2*R*F < W <= 4*R*F` 区间过早启用半网格剪枝。

当一对仍处于搜索范围内的颗粒落在被剪掉的相邻单元中时，接触候选会被
漏掉。串行和 DMP 的粒子所有权、ghost 布局及访问顺序不同，因此同一输入
还可能表现为 rank-dependent 接触集合、力和轨迹差异。

### 2.2 本项目如何处理

`code/desgrid_mod.f:1116-1126` 在 `USE_JKR=.TRUE.` 时把阈值系数由 2 改为
4；`USE_JKR=.FALSE.` 时保留原生系数 2，以维持明确要求的原 MFIX 兼容
基线。也就是说，当前改动修复了本项目 JKR 路径，不声称已经普遍修复原生
Hertz/LSD。

固定 8-rank、2×2×2 双尺度工况对此作了定量验证：`W=3.333 um` 位于错误
区间时，两档 JKR 串并行一致，而 JKR-off 负对照稳定出现 serial 2 对、DMP
4 对；`W=5.0 um > 4.8 um` 时 JKR 与 JKR-off 均一致。源码门控契约位于
`validation/tests/test_next_stage_contract.py:182-199`，数值证据详见
`validation/docs/VALIDATION.md` 的“固定 8-rank 双尺度颗粒群”。

## 3. DMP 的本地 CC 不能单独承担跨 rank 成对历史身份

### 3.1 Global ID 不会因跨分区而切换

MFIX 支持论坛对 DMP 粒子身份的说明是：本地索引 `LL` 会变化，
`iGlobal_ID(LL)` 不会变化；global particle ID 会一直附着于同一颗粒，直到
颗粒离开计算域。参见 [Particle ID numbers][particle-id]。因此“颗粒跨 MPI
分区时 Global ID 切换”不是本次角速度差的原因。需要区分三件事：

- 同一个输入文件分别启动 serial 和 DMP 时，初始化/分区过程可能给物理颗粒
  不同的初始 global ID；这属于两次运行之间的标号差异。
- 同一次 DMP 运行内，颗粒迁移后 rank、本地数组位置 `LL` 和接触槽 `CC` 可以
  改变，但该颗粒的 `iGlobal_ID` 保持不变。
- 有颗粒入口/出口时，新粒子的 ID 分配和离域粒子的 ID 生命周期另有规则；
  这不等于域内迁移时给现存颗粒改号。

MFIX 26.1 用户手册只把 `VTK_PART_ID` 描述为输出 particle ID，并给出 DMP
运行方法；目前未找到它对 `iGlobal_ID` 跨分区不变性或 `CC/PFT_NEIGHBOR`
迁移语义的用户级保证。上述身份结论来自官方支持答复和 26.1.2 源码，而不是
手册中一个已发布的接触历史契约。

### 3.2 原生 PFT 不会遇到本次 PFR 的编号符号错误，但有独立连续性风险

原生 `model/des/calc_force_dem.f:385-408` 直接用
`PFT_NEIGHBOR(:,CC)` 积分 `VREL_T`，没有按 global-ID 大小给物理增量乘正负号，
所以不会复现本项目旧 PFR 的“交换编号后滚动力矩反转”。PFT 是
owner-directed history；若 owner/neighbor 对调，`VREL_T` 和 PFT 都应反号。

但原生 26.1.2 的跨 rank 连续性是另一问题：
`model/des/neighbour.f:77-104` 只按同一 owner 下的本地 neighbor index 恢复
PFT，找不到就清零；`model/des/mpi_pack_des_mod.f:485-520` 把 pair 循环设成
零次，并注释了 PFT 打包调用。MFIX 支持论坛 2025-12 至 2026-01 的
[DMP contact-time 讨论][contact-time]也复现了存入 `PFT_NEIGHBOR(:,CC)` 的
成对接触时间在分区边界重置；开发者确认 DMP 的 CC 逻辑存在问题，并建议
维护可通信的 pair identity，而不是跨 rank 直接依赖 CC。

该证据说明：对任何需要严格连续的成对历史，原生 rank-local CC 本身都不是
充分的跨分区身份。它不等价于证明所有原生 Hertz/LSD 工况都会出现可见的
切向力误差；本项目目前没有把“未经修改的 26.1.2 PFT 跨 owner 定量失效”做成
独立负对照。因此这里记录为源码可见且获官方支持确认的上游限制/风险，并保留
这一证据边界。

### 3.3 本项目处理和本次缺陷归属

本项目没有另建平行接触表，而是在 MFIX 的 neighbor CC 生命周期上使用
global-ID pair 搬运历史：DMP/restart 先按两个 global ID 暂存，官方邻居重建
产生新 CC 后恢复；owner 反转时 owner-directed PFT 反号，pair-invariant PFR
不反号。active migration 和 DMP restart 已用双精度事件验证这条链路。

本次逐分量角速度差应归类为本项目旧滚动实现的物理/参考系错误：旧代码用
global-ID 大小给 PFR 增量定向，同时直接把切平面位移当力矩轴，导致相同物理
初态仅重标号就改变力矩。修复是取消编号参与物理更新，统一 PFR 的
pair-invariant 定义和 `normal x PFR` 力矩映射。它应写入迁移和验证缺陷记录，
不能列为 MFIX 上游 limitation。这里保留本节，是为了说明为什么本项目仍需
扩展原生 PFT 的 DMP history 搬运机制。

## 4. 状态与后续决策

前两项不足已有本项目可重复负对照；第三项有 26.1.2 源码证据、官方支持帖
复现和开发者答复，本项目的 JKR history 路径另有正向迁移/restart 验证，但
尚未增加纯原生 PFT owner-migration 负对照。所有 JKR 修复均保持 MFIX 原有
槽或邻居搜索结构。若未来要把修复推广到 JKR-off，应先补齐对应负对照并改为
物理连续性验收，再单独评估改变 MFIX 原生结果的兼容性影响。

[particle-id]: https://mfix.netl.doe.gov/forum/t/particle-id-numbers/3299
[contact-time]: https://mfix.netl.doe.gov/forum/t/abnormal-particle-contact-time-in-dem-simulations-with-dmp-parallelization/6581
