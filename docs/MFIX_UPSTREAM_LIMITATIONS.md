# MFIX 原生实现的已知不足与本项目处理

更新：2026-08-24

本文记录迁移和验证过程中由源码审查及定量工况确认的 MFIX 原生行为。
这里的“原生”指本仓库 `model` 中的 MFIX-26.1.2 基线；`code` 是本项目的
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

## 3. 状态与后续决策

两项不足都已有可重复的负对照，并且 JKR 修复均保持 MFIX 原有槽或邻居
搜索结构。若未来要把修复推广到 JKR-off，应先把负对照改为物理连续性或
串并行一致性验收，再单独评估改变 MFIX 原生结果的兼容性影响。
