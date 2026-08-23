# JKR Validation Test Progress

更新：2026-08-23

## 当前基线

- 分支：`codex/jkr-validation`。
- MFIX：26.1.2。
- 本阶段父 commit：`a6f4a803f5eafd102c0ffd70c1f5064d39dc3054`。
- 串行 solver SHA-256：`1cf3d9f4731a58174320db718b014bcf84a03c6b1e9ed3e287eea5a59992eb76`。
- DMP solver SHA-256：`003290b15d19c1d1edee39c5063181935c024dc8c0b8f2b55d62a65098e03d81`。
- 当前结果快照：`/home/yuxing/mfix/MfixVer26/validation/results`；每次统一运行替换该快照，`RUN_INFO.txt` 记录时间与版本。
- 已通过：91 个单元测试、29 个 canonical case 字节同步、15 个 serial JKR 工况和 5 个固定 2-rank DMP 工况。

## 本轮目标与顺序

1. 已完成：29 个 Git 跟踪 canonical case、显式 e=1/e=0.6 档案及逐字节同步检查。
2. 已完成：e=1 固定 0.95/1.05 JKR、对应 Hertz、独立 RK4 参考、CSV 和四类 contact 图。
3. 已完成：异质表面能的均质退化、不同粒子 gamma、ID 交换、low-high 对称、wall 和默认回退。
4. 已完成：2-rank DMP 跨区新接触、活跃历史迁移、wall、cloud，按 global ID/pair ID 比较。
5. 已完成接口：默认关闭的双精度 JKR history 诊断；正式 restart 边界矩阵仍为下一阶段。
6. 较后：壁面 single face/edge/corner/facet seam；多壁面同时接触不作为近期阻塞项。

## 结果目录组织约定

后续新运行参考 MFIX23 的可浏览方式，按
`测试类别 / 物理子类 / case / 执行模式` 组织。例如：

- `contact_physics/pair_normal/pair_normal_capture/serial`；
- `contact_physics/pair_normal/pair_gamma_low_high/serial`；
- `dmp/pair_cross/pair_cross_jkr/{serial,dmp}`；
- `restart/pair_history/restart_pair_pfr/{continuous,restart1,restart2}`。

最外层固定为 `validation/results`。每次 runner 调用先清理上一份结果快照，
再写入一套完整结果；需要复核旧结果时，从相应 Git commit、canonical case
和已记录的 solver hash 重新运行。`dmp` 表示 DMP 执行方式，进程数不写入
目录名，而在 `RUN_INFO.txt` 和 metadata 中显式记录；下一阶段固定为
`mpi_processes=2`。完整契约见 `docs/NEXT_STAGE_PROMPT.md`。

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
| 验证框架与统一 restart | PASS | 91 unit tests；restart 证据由统一 all 运行重建 |
| canonical case 资产 | PASS | 29 case、58 个输入文件与集中生成器逐字节一致 |
| 文献和旧测试过程审查 | PASS | 4 篇直接相关论文复制到 `docs/references`，旧目录只读 |
| 原生 Hertz 保护 | PASS | 3 个 e=0.6 原生保护和 4 个 e=1 对应速度 case；官方/迁移轨迹及双精度末态一致 |
| JKR 理论关键点 | PASS | 已单测最大拉力、临界分离、零重叠和平衡点，明确最大拉力与临界分离不是同一点 |
| 基础 capture/rebound | PASS | e=1 wall/pair：0.95× 捕获、1.05× 断颈回弹；e=0.6 case 明确改名；逐点力、CSV 和四类图齐全 |
| 速度扫描 | DEFERRED | 按当前要求不搜索转变区间；后续若需要只增加显式固定速度矩阵 |
| 异质表面能 | PASS | 7 个 e=1、0.05R case；几何平均、交换/ID 对称、wall 与 nonpositive fallback |
| 2-rank DMP | PASS | pair cross/new/active migration、wall cross、cloud；2 ranks、2×1×1 |
| restart 内部历史 | INTERFACE READY | `MFIX_JKR_DIAGNOSTICS` 默认关闭；正式 restart 边界矩阵待执行 |
| 壁面几何 | DEFERRED | 多壁面同时接触后移；先不阻塞异质和 DMP |

## 已知范围

- 普通 Hertz/JKR approach case 的初始表面间距统一为 `0.05R = 5e-8 m`；wall 粒心 y 为 1.05 μm，pair 粒心 y 为 8.975/11.025 μm。
- restart case 为建立 PFT/PFR/state 历史，保留其专用初始接触或正重叠，不套用 0.05R。
- 预期单颗粒不会同时接触超过 8 个壁面；超过 8 个墙面接触暂不纳入。
- SuperDEM/GSP 不在当前 JKR 验证范围。
- 标准 VTP 会把双精度值转换为 Float32；restart 内部历史必须使用独立、默认关闭的双精度诊断输出。
- 当前 wall contact 历史槽初值为 8；JKR+DMP 超过容量会显式报错，本阶段不覆盖该范围。
- DMP unpack 的 neighbor fallback 分支已只读审查；五个真实 2-rank case 未产生定位失败证据，因此未修改该分支。
- active migration 的严格首差出现在 ownership 迁移时刻；双精度诊断证明 PFT/PFR/JKR state 的 queue/restore 连续。
- 诊断不复用 `DES_USR_VAR(2:)`，无固定 case、时间或粒子 ID；关闭时不创建输出。
- 大型运行输出不提交 Git，只记录绝对路径、时间和 solver hash。
