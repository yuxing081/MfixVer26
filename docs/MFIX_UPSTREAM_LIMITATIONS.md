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

## 3. 颗粒跨 DMP 分区时会丢失 PFT 切向历史

### 3.1 错误发生的历史过程

原生 MFIX 把颗粒对的切向位移保存在 `PFT_NEIGHBOR(:,CC)` 中。其中 `CC` 是
当前 rank 邻居表里的本地接触槽，并不是跨 rank 稳定的物理接触身份。当仍在
接触的一对颗粒跨越 MPI 分区，或迁移后接触 owner 反转时，会经历以下过程：

1. 迁移前，PFT 只存在于旧 rank 的旧 `CC` 槽中。
2. 颗粒迁移后，其 rank、本地粒子索引 `LL` 和接触槽 `CC` 可以改变；接触对
   也可能由 `(i,j)` 改为 `(j,i)`。颗粒 Global ID 本身不会因此改变。
3. 原生 26.1.2 的 particle pair packet 没有形成有效搬运闭环：
   `model/des/mpi_pack_des_mod.f:485-520` 令 pair 循环为零次，并注释掉 PFT
   的 pack 调用。
4. 随后的邻居表重建只在同一旧 owner 的本地 neighbor index 中查找历史，见
   `model/des/neighbour.f:77-104`。如果新 `CC` 无法匹配旧槽，PFT 被初始化为
   零，连续接触的切向弹簧历史因此重新开始。

这不是 Global ID 在跨区时切换。MFIX Support 的
[Particle ID numbers][particle-id] 明确说明，DMP 中改变的是 `LL`，
`iGlobal_ID(LL)` 保持不变并随颗粒移动。MFIX 26.1 用户手册只说明 particle ID
输出和 DMP 运行方式，目前未给出 `CC/PFT_NEIGHBOR` 跨分区连续性的用户级
保证。

MFIX Support 的 [Abnormal particle contact time in DEM simulations with
DMP parallelization][contact-time] 提供了直接佐证：提问者把成对接触时间存入
`PFT_NEIGHBOR(4,CC)`，观察到它在 DMP 分区边界重置，而串行和单 rank DMP
正常；开发者确认 DMP 的 `CC` 逻辑存在问题，并建议维护可跨 rank 传递的接触
pair identity，不要直接依赖 `CC`。该案例保存的是接触时间而不是原生三分量
PFT，但触发的是相同的 `CC` 身份和迁移缺口；原生 PFT pack 被禁用则由上述
26.1.2 源码直接确认。

### 3.2 本项目如何处理 PFT

本项目复用 MFIX 已有的 neighbor、particle migration 和 pair packet 机制，
没有引入第二套长期接触表：

1. 迁移时把旧 owner 和 neighbor 的两个 Global ID 作为稳定的物理 pair
   identity，并随 pair packet 发送 PFT。
2. 接收端不在 particle unpack 阶段猜测新的 `CC`，而是先按 Global-ID pair
   放入 pending queue。
3. 等 MFIX 完成邻居表重建、产生有效的新 `CC` 后，再将历史恢复到该槽。
4. 若 owner/neighbor 顺序反转，PFT 同时取负，因为 PFT 与相对切向速度一样，
   是按 `owner -> neighbor` 定义的有方向历史。
5. 排序、普通邻居重建和 restart 也使用相同的 pair identity 与方向规则，
   避免不同入口各自解释接触历史。

23 版 JKR 修改已经启用了 pair PFT 的 pack/unpack，并通过 `ADD_PAIR` 恢复到
接收 rank 的接触槽；26 版迁移保留了这一意图，但适配了 MFIX 26.1.2 新的邻居
重建和 restart 生命周期，改为 Global-ID pending restore。非零 PFT 的 active
owner migration 和 DMP restart 已由双精度事件验证。

### 3.3 PFR 和 JKR 接触状态的相同点、不同点与处理

新增 PFR 滚动历史和 JKR tensile-neck 接触状态与 PFT 的相同点是：三者都
属于同一个物理颗粒对，必须共享同一 Global-ID pair、同一权威 `CC`、同一
pack/unpack、pending restore、邻居重建和 restart 生命周期。接触真正断开时，
三者也在同一接触槽中一起清除，不能因 rank 或 owner 改变而提前重置。

不同点在 owner 反转时的变换规则：

| 历史量 | owner 反转时 | 原因 |
| --- | --- | --- |
| PFT | 取反 | owner-directed 切向位移，相对切向方向随 owner 反转 |
| PFR | 保持不变 | Chen 滚动位移定义为 pair-invariant 接触平面向量；法向反转后力矩自然反向 |
| JKR contact state | 保持不变 | tensile neck 是否激活是整个物理接触对的无方向逻辑状态 |

因此三者一起搬运、一起按 Global-ID pair 找回新 `CC`，但只对 PFT 做反号。
本次发现的逐分量角速度差正是旧 PFR 实现错误地使用了编号定向并在 owner
反转时变号；它属于新增滚动模型的物理/参考系缺陷，不是原生 PFT 跨区丢失
本身。现在 PFR 不再依赖 Global-ID 大小，并通过 `normal x PFR` 映射滚动力矩。
墙面 PFT/PFR/state 则作为粒子附着历史随 particle packet 搬运，继续使用
MFIX 原有 facet/history slot。

## 4. 状态与后续决策

前两项不足已有本项目可重复负对照；第三项有 26.1.2 源码证据、MFIX Support
的同类 `PFT_NEIGHBOR(:,CC)` 历史跨区重置复现和开发者答复，本项目的 JKR
history 路径另有正向迁移/restart 验证，但尚未增加纯原生三分量 PFT 的
owner-migration 负对照。所有 JKR 修复均保持 MFIX 原有
槽或邻居搜索结构。若未来要把修复推广到 JKR-off，应先补齐对应负对照并改为
物理连续性验收，再单独评估改变 MFIX 原生结果的兼容性影响。

[particle-id]: https://mfix.netl.doe.gov/forum/t/particle-id-numbers/3299
[contact-time]: https://mfix.netl.doe.gov/forum/t/abnormal-particle-contact-time-in-dem-simulations-with-dmp-parallelization/6581
