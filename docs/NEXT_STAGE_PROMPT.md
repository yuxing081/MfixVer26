# MFIX 26.1.2 JKR 当前验证阶段

更新日期：2026-08-24。

旧的异质表面能、固定 2-rank DMP 和 restart 边界实施计划已完成并归档到
`docs/archive/NEXT_STAGE_PROMPT_2026-08-23.md`。开始工作前检查生产仓库与
嵌套 `validation/` 仓库的实际 Git 状态，不得重置或覆盖用户修改。

## 当前基线

- 球形颗粒 JKR、PFT/PFR、逐颗粒表面能、邻居/DMP 历史和 DES RES 2.0
  已迁移；SuperDEM/GSP 与超过默认 8 个 wall history slots 不在范围内。
- serial restart 和固定 2 ranks、2×1×1 DMP restart 已通过；诊断由
  `MFIX_JKR_DIAGNOSTICS` 显式启用，默认关闭，本阶段暂不移动。
- 法向 capture/rebound 已与独立 JKR 理论比较，不重复增加 tensile-neck
  法向测试。

## 本阶段顺序

1. 已完成：定量验证颗粒对和颗粒–墙面的斜碰、滑动与滚动；独立公式、
   PFT/PFR、切向力、重建滚动力矩、速度/角速度和 stick/slip 覆盖均通过。
2. 已完成：使用真实 ASCII STL 验证 wall single face、edge、corner 和
   triangulated seam；接触计数、合力、facet identity 与初始历史均通过，
   并修复 JKR 共面 seam 双计。
3. 已完成：固定 2×2×2、8-rank 双尺度颗粒群。JKR 两尺度和高尺度
   JKR-off 均与单核逐字段一致；低尺度 JKR-off 稳定复现原始系数 2 的
   rank-dependent 漏接触负对照，直接证明 JKR 系数 4 修正有效。
4. 已完成：`./run_validation.sh all` 共 42/42 PASS；117 项 unit 中
   113 PASS、4 个可选 `fparser` 测试跳过。两个仓库的结果来源 commit
   均为 clean。

本文件列出的本阶段工作已全部完成。后续仍明确不在当前验收范围的内容为：
更多 rank/空间分解遍历、超过默认 8 个 wall history slots、SuperDEM/GSP、
能量分解，以及自动临界速度搜索；除非用户另行扩展范围，不将其恢复为
当前待办。

每项必须先定义固定参数和验收判据，再增加最小 case、独立比较和真实运行。
完成一项后立即更新 `docs/MIGRATION.md`、`validation/docs/VALIDATION.md`、
`validation/docs/TEST_PROGRESS.md`，分别提交受影响的 Git 仓库。

## DMP 粗网格专项解释

令 DES search range 为当前邻居搜索使用的最大距离尺度。低尺度 case 的
cell width 必须小于其两倍，但仍落在 MFIX 官方系数 2 与 JKR 修正系数 4
会产生不同 half-cell 判定的区间；高尺度 case 的 cell width 大于其两倍。
两个 case 使用可重复的颗粒群和显式跨 rank 接触，不用随机碰撞是否发生来
决定测试有效性。

验收时 JKR serial 只与同参数 JKR DMP 比较，JKR-off serial 只与同参数
JKR-off DMP 比较。粒子按 global ID、接触按物理 pair ID 对齐，记录轨迹、
接触集合、PFT/PFR/state、首次差异和最大误差；不得把 JKR 与 Hertz 的物理
轨迹相互比较作为一致性要求。

实际固定参数与结果：`W=3.333/5.0 μm`，DES search range 为 `2.4 μm`，
两倍阈值为 `4.8 μm`。每个 case 初始 owner 覆盖 rank 0–7，并形成四组
跨 x pair。低尺度 JKR-off 被定义为预期负对照：首次接触集合差异为
serial 2 对、DMP 4 对；其余三项 serial–DMP 的 x/v/omega/force 最大误差
均为零。该纯法向专项中 PFT/PFR 预期为零；非零历史由既有迁移/restart
工况覆盖。
