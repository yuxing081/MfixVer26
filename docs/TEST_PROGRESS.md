# JKR Validation Test Progress

更新：2026-08-22

## 当前基线

- 分支：`codex/jkr-validation`。
- MFIX：26.1.2。
- 生产源码基线：`954ea0b9639592b9632c4f8fba92d23eabf42b99`。
- 串行 solver SHA-256：`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`。
- DMP solver SHA-256：`3523a0d4352431534f8576e905201c04abc7105ba6902f721552916a182ab26c`。
- 最后完整证据：`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T095448.025415Z-107484`。
- 已通过：40 个单元测试、serial pair/wall RESTART_1 末态逐字节回归、RESTART_2 敏感性负对照。

## 本轮目标与顺序

1. 原生 Hertz 保护：官方 MFIX26 与 JKR solver 在 `use_jkr=.False.` 下对照。
2. JKR 理论关键点与基础捕获/回弹：平衡点、零重叠稳定点、最大拉力/拉脱点、临界分离，随后运行 wall/pair capture 与 rebound。
3. wall/pair 速度扫描：先粗扫定位转变区间，再加密一次，输出统一 CSV 和五类图。
4. 异质表面能：均质退化、不同粒子 gamma、ID 交换、low-high 对称、wall、默认回退。
5. 壁面几何：single face、edge、corner、facet seam。
6. 2-rank DMP：cross-rank 新接触、活跃历史迁移、wall、cloud，并按 global ID/pair ID 比较。
7. restart 内部边界诊断：NEIGHBORS、PFT、PFR、JKR state、facet 和 wall histories 的双精度逐字段连续性。

## 当前判据

- 不以 solver 退出码 0 作为通过。
- 理论测试必须报告关键点误差和符号约定。
- 原生保护必须与未修改 MFIX26 solver 比较。
- DMP 必须按 global particle ID 和 canonical contact pair ID 对齐。
- restart 必须比较 continuous、RESTART_1、RESTART_2；内部历史使用双精度。
- MFIX23 结果只作过程参考，不作为真值。

## 活动记录

| 阶段 | 状态 | 最近证据/阻塞 |
|---|---|---|
| 验证框架与统一 restart | PASS | 40 unit tests；最终 all 证据见当前基线 |
| 文献和旧测试过程审查 | IN PROGRESS | 只读审查 `versions/docs/references` 与 `versions/test` |
| 原生 Hertz 保护 | NEXT | 需要定位/构建未修改 MFIX26 基准 solver |
| JKR 理论关键点 | NEXT | 独立参考函数已有基础，待补关键点单测 |
| 基础 capture/rebound | PENDING | 等关键点和 Hertz 基线 |
| 速度扫描 | PENDING | 等基础工况 |
| 异质表面能 | PENDING | 等基础工况 |
| 壁面几何 | PENDING | 等基础工况 |
| 2-rank DMP | PENDING | serial 通过后开始 |
| restart 内部历史 | PENDING | 需要最小可选双精度诊断接口 |

## 已知范围

- 预期单颗粒不会同时接触超过 8 个壁面；超过 8 个墙面接触暂不纳入。
- SuperDEM/GSP 不在当前 JKR 验证范围。
- 大型运行输出不提交 Git，只记录绝对路径、时间和 solver hash。
