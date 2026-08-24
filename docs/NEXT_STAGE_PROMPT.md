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

1. 定量验证颗粒对和颗粒–墙面的斜碰、滑动与滚动；比较 PFT/PFR、切向力、
   滚动力矩、速度和角速度，不只检查趋势或 solver 退出码。
2. 验证 wall single face、edge、corner 和 triangulated seam；检查接触计数、
   合力、facet identity，以及 PFT/PFR/state 在 facet 切换时的连续性。
3. 建立 2×2×2、8-rank 颗粒群工况，与单核运行比较。使用两个计算域/网格
   尺度：cell width 分别位于 `2*DES search range` 阈值两侧；JKR 和
   JKR-off 原生 MFIX 分别进行 serial–DMP 对照。
4. 前三项完成后再执行完整综合回归和阶段总结。

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
