# JKR Validation Test Progress

更新：2026-08-23

## 当前基线

- 分支：`codex/jkr-validation`。
- MFIX：26.1.2。
- 生产源码基线：`954ea0b9639592b9632c4f8fba92d23eabf42b99`。
- 串行 solver SHA-256：`0e7a4e7e8fad556481bcb11bd76aed0508b180c7a5e0944aaca681dfc8a0ab6f`。
- DMP solver SHA-256：`3523a0d4352431534f8576e905201c04abc7105ba6902f721552916a182ab26c`。
- 最新统一全套证据：`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T151407.503647Z-5944`（unit、7 Hertz、8 JKR normal、2 restart 共 18 项全部 PASS）。
- 前一 restart 基线：`/home/yuxing/mfix/MfixVer26/validation/results/runs/20260822T095448.025415Z-107484`。
- 已通过：67 个单元测试、17 个 canonical case 字节同步、7 个官方 MFIX26 Hertz 对照、8 个 serial JKR normal 工况、32 张 contact 图，以及既有 serial restart 回归。

## 本轮目标与顺序

1. 已完成：17 个 Git 跟踪 canonical case、显式 e=1/e=0.6 档案及逐字节同步检查。
2. 已完成：e=1 固定 0.95/1.05 JKR、对应 Hertz、独立 RK4 参考、CSV 和四类 contact 图。
3. 下一阶段：异质表面能的均质退化、不同粒子 gamma、ID 交换、low-high 对称、wall 和默认回退。
4. 随后：2-rank DMP 跨区新接触、活跃历史迁移、wall、cloud，按 global ID/pair ID 比较。
5. 随后：restart 边界 NEIGHBORS、PFT、PFR、JKR state、facet 和 wall histories 双精度诊断。
6. 较后：壁面 single face/edge/corner/facet seam；多壁面同时接触不作为近期阻塞项。

## 结果目录组织约定

后续新运行参考 MFIX23 的可浏览方式，按
`测试类别 / 物理子类 / case / 执行模式` 组织。例如：

- `contact_physics/pair_normal/pair_normal_capture/serial`；
- `contact_physics/pair_normal/pair_gamma_low_high/serial`；
- `dmp/pair_cross/pair_cross_jkr/{serial,dmp2}`；
- `restart/pair_history/restart_pair_pfr/{continuous,restart1,restart2}`。

最外层仍为不可覆盖的 `validation/results/runs/<timestamp>`。既有历史结果
不移动、不改名；新布局只应用于后续运行。完整契约见
`docs/NEXT_STAGE_PROMPT.md`。

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
| 验证框架与统一 restart | PASS | 67 unit tests；restart 证据见当前基线 |
| canonical case 资产 | PASS | 17 case、34 个输入文件与集中生成器逐字节一致；15 个普通 approach case 均为 0.05R 初始间距 |
| 文献和旧测试过程审查 | PASS | 4 篇直接相关论文复制到 `docs/references`，旧目录只读 |
| 原生 Hertz 保护 | PASS | 3 个 e=0.6 原生保护和 4 个 e=1 对应速度 case；官方/迁移轨迹及双精度末态一致 |
| JKR 理论关键点 | PASS | 已单测最大拉力、临界分离、零重叠和平衡点，明确最大拉力与临界分离不是同一点 |
| 基础 capture/rebound | PASS | e=1 wall/pair：0.95× 捕获、1.05× 断颈回弹；e=0.6 case 明确改名；逐点力、CSV 和四类图齐全 |
| 速度扫描 | DEFERRED | 按当前要求不搜索转变区间；后续若需要只增加显式固定速度矩阵 |
| 异质表面能 | NEXT | serial 均质退化、异质组合、ID 对称和 fallback |
| 2-rank DMP | PENDING | 异质表面能后，先 pair cross/new/migration |
| restart 内部历史 | PENDING | DMP 后增加最小可选双精度诊断接口 |
| 壁面几何 | DEFERRED | 多壁面同时接触后移；先不阻塞异质和 DMP |

## 已知范围

- 普通 Hertz/JKR approach case 的初始表面间距统一为 `0.05R = 5e-8 m`；wall 粒心 y 为 1.05 μm，pair 粒心 y 为 8.975/11.025 μm。
- restart case 为建立 PFT/PFR/state 历史，保留其专用初始接触或正重叠，不套用 0.05R。
- 预期单颗粒不会同时接触超过 8 个壁面；超过 8 个墙面接触暂不纳入。
- SuperDEM/GSP 不在当前 JKR 验证范围。
- 标准 VTP 会把双精度值转换为 Float32；restart 内部历史必须使用独立、默认关闭的双精度诊断输出。
- 当前 wall contact 历史槽初值为 8；JKR+DMP 超过容量会显式报错，本阶段不覆盖该范围。
- DMP unpack 的邻居 fallback 分支使用 owner IJK 定位 neighbor，已列为 2-rank 跨区测试重点观察项，尚不先行修改生产代码。
- 大型运行输出不提交 Git，只记录绝对路径、时间和 solver hash。
