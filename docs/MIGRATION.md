# JKR 向 MFIX 26.1.2 迁移记录

## 当前状态

- 阶段：第五阶段固定 2-rank DMP restart 已验收；历史章节保留各阶段当时结论。
- 模型范围：仅球形颗粒；SuperDEM/GSP 异形颗粒不做 JKR 适配。
- restart 原则：采用 MFIX 26.1.2 官方修复和 global-ID 映射，PFR/JKR state/wall history 仿照官方 PFT 生命周期扩展，不自行开发另一套邻居 restart 映射。
- 验证边界：能量统计、DES_USR_VAR(2:31) 和历史 seed 只放在 validation 的具体 case 或 analysis 工具中。
- 迁移开始日期：2026-08-22。
- 新项目根目录：/home/yuxing/mfix/MfixVer26。
- MFIX 26.1.2 官方安装源码：/home/yuxing/miniforge3/envs/mfix-26.1.2/share/mfix/src。
- 优先基准快照：/home/yuxing/mfix/MfixVer26/model；与安装源码的 model/ 逐文件比较一致。
- 旧 MFIX 23.4.1 基准：/home/yuxing/miniforge3/envs/mfix-23.4.1/share/mfix/src。
- 旧 JKR 仓库：/home/yuxing/mfix/P3/JKR/versions/code。
- 旧 JKR 分支/提交：jkr-cc-state-refactor，ea55027498a4b939a2982cbb77d137201b034ffc（Expand the useful range of γ）。
- 旧工作树：dirty；仅 calc_force_dem.f 有未提交修改（18 additions、3 deletions），内容为 restart 回归专用 PFR seed 和 DES_USR_VAR(23:31) 力/力矩诊断，不属于生产物理。
- 已核验工具链：MFIX 26.1.2、GNU Fortran 15.3.0、CMake 4.2.3；构建入口为激活 mfix-26.1.2 环境后的 build_mfixsolver。mpifort wrapper 存在，但必须在完整激活的 conda 环境中再次验证 DMP 构建。

本项目只维护相对于 MFIX 26.1.2 修改或新增的文件。不得把 23.4.1 文件整份覆盖到 26.1.2，也不得修改 conda 环境中的官方源码。

## 第一阶段结论

26.1.2 的接口变化较大，不能直接开始移植：颗粒和墙面接触力路径已加入 GSP 等新版逻辑；restart 的邻居身份由 23.4.1 的二维接触对改为一维 owner-neighbor 全局 ID 映射；DMP pair packet 虽保留框架，但 pair 循环、ADD_PAIR 和 PFT pack/unpack 在官方 26.1.2 中仍被禁用或注释。第二阶段只能复制下表中确定需要修改的 26.1.2 文件，第三阶段必须逐模块语义重写和编译验证。

## 文件迁移矩阵

| 旧文件或旧意图 | 26.1.2 对应文件 | 修改目的 | 上游状态 | 迁移决定 | 风险 |
| --- | --- | --- | --- | --- | --- |
| allocate_arrays.f | model/allocate_arrays.f | 释放 PFR、pair/wall JKR state 和 wall PFR | 新版仅处理原生数组 | 在新版 DEM deallocate 路径增量追加 | 中：生命周期遗漏 |
| calc_collision_wall_mod.f | model/des/calc_collision_wall_mod.f | 墙面 JKR、负重叠颈部、PFT/PFR、滚动力矩、facet state | 新版有 STL、CFRELVEL_WALL、GSP、Hertz/VDW，无 JKR/PFR/state | 基于新版接口语义重写；诊断拆出 | 极高：facet slot、投影和 GSP 分支 |
| calc_force_dem.f | model/des/calc_force_dem.f | Chen/JKR 法向力、异质表面能、PFT/PFR、滚动力矩、neck state | 新版已重构 Hertz/VDW/GSP，无 JKR/PFR/state | 完全语义重写；USE_JKR=.FALSE. 必须走未改动的原生 Hertz/LSD 路径 | 极高：物理分支和新版数据流 |
| cfnewvalues.f | model/des/cfnewvalues.f | 旧版仅调用 energy audit step | 新版为 GSP state/integrator | 生产代码采用官方文件；审计调用放 validation-only 构建 | 低 |
| check_solids_dem.f | model/check_data/check_solids_dem.f | JKR 参数、组合和等效量检查 | 新版只有原生 collision/cohesion/Hamaker 检查 | 在新版检查框架增量重写 | 中 |
| des_allocate_mod.f | model/des/des_allocate_mod.f | 分配/grow PFR、pair/wall state；ADD_PAIR 避免重复 slot | 新版仅分配 PFT、wall facet/PFT 和 SuperDEM 历史 | 扩展新版 allocate/grow/deallocate/ADD_PAIR | 高：CC 重排和重复接触 |
| des_energy_audit.f | 无官方对应 | 接触势能及阻尼、滑移、滚动损失的 MPI 审计 | 上游无实现 | 不进入生产 code/；移至 validation/analysis，仅显式验证构建启用 | 中 |
| des_init_namelist.f | model/des/des_init_namelist.f | USE_JKR、表面能、临界滚动角、安全系数及默认值 | 新版仅有 Hamaker/VDW 参数 | 在新版默认值路径追加；默认 USE_JKR=.FALSE. | 中 |
| desgrid_mod.f | model/des/desgrid_mod.f | half-cell 搜索阈值和跨 rank 对称接触对 | 旧 4*MAX_RADIUS 修正未上游化；26.1.2 仍为 2*MAX_RADIUS，单侧过滤仍存在 | 在新版邻居构建上语义重做并设专项 DMP 测试 | 极高 |
| desnamelist.inc | model/include/desnamelist.inc | 声明 JKR 输入参数 | 26.1.2 无 JKR 声明 | 仅追加所需声明，不复制旧 include | 低 |
| discretelement_mod.f | model/des/discretelement_mod.f | JKR 派生矩阵、逐颗粒表面能访问、PFR、pair/wall state | 新版仅有 PFT、wall facet/PFT、Hamaker/phase surface energy | 增量加入明确的数据布局和访问函数 | 高 |
| mpi_init_des_mod.f | model/des/mpi_init_des_mod.f | particle/pair packet 尺寸 | 新版 pair size 为 11，SuperDEM 另加 8；无 JKR/wall history | 依据唯一字段清单重新计算 packet size | 极高 |
| mpi_pack_des_mod.f | model/des/mpi_pack_des_mod.f | pair PFT/PFR/state 与粒子附着 wall history 打包 | 新版 pair 循环实际停用，PFT pack 被注释 | 保留新版迁移框架，重新启用 pair 协议并增量打包 | 极高 |
| mpi_unpack_des_mod.f | model/des/mpi_unpack_des_mod.f | 对称恢复 pair/wall history 和 ADD_PAIR slot | 新版 ADD_PAIR/PFT unpack 被注释 | 与 pack 使用同一协议表语义重写 | 极高 |
| neighbour.f | model/des/neighbour.f | 邻居重建时恢复 PFT/PFR/state；owner 反转时向量反号 | 新版仅同向恢复 PFT 和 SuperDEM history | 采用新版重建，再追加 PFR/state 及 reversed-owner 规则 | 极高 |
| 无旧独立文件；由历史迁移需求发现 | model/des/des_functions.f | 粒子排序/空间重排时 remap CC 历史 | 新版只 remap PFT 和 SuperDEM history | 新增到需修改清单；扩展所有排序/exchange remap 路径 | 极高 |
| read_namelist.f | model/read_namelist.f | 旧差异只有末尾空白/换行 | 无需上游替代 | 废弃旧修改；默认不复制，只有编译证明 parser 需要挂接时才最小修改 | 低 |
| read_res0_des.f | model/des/read_res0_des.f | 读取 PFR、pair/wall state、facet 和 wall PFT/PFR | 新版已读取 global IDs、neighbor/PFT，并由官方 map 重建邻居 | 在新版布局和 map 上追加版本化字段 | 极高 |
| write_res0_des.f | model/des/write_res0_des.f | 写入 PFR、pair/wall state、facet 和 wall PFT/PFR | 新版写 global neighbor ID/PFT；不写旧式 NEIGHBOR_INDEX | 在新版布局末尾追加字段，不迁移旧写法 | 极高 |
| 23.4.1 read_res1_des_mod.f restart 修补意图 | model/des/read_res1_des_mod.f | global-to-local collision owner/neighbor map | 26.1.2 已实现 READ_PAR_COL → MAP_cARRAY_TO_PROC → GLOBAL_TO_LOC_COL；23.4.1 对应循环被 FIXME 注释 | 完全采用官方 26.1.2；不复制到 code/ | 低 |
| 旧式二维 iPAR_COL 和直接保存 NEIGHBOR_INDEX | 新版一维 neighbor global ID + 重建 index | 旧 restart 接触身份格式 | 已被新版官方设计替代 | 废弃 | 高 |

## 第二阶段复制到 code/ 的文件

以下 16 个 MFIX 26.1.2 文件已按 basename 复制到扁平 code/ 顶层。来源路径保留如下；扁平布局是 MFIX CMake 的 UDF source-replacement 规则要求：

1. model/allocate_arrays.f
2. model/check_data/check_solids_dem.f
3. model/include/desnamelist.inc
4. model/des/calc_collision_wall_mod.f
5. model/des/calc_force_dem.f
6. model/des/des_allocate_mod.f
7. model/des/des_functions.f
8. model/des/des_init_namelist.f
9. model/des/desgrid_mod.f
10. model/des/discretelement_mod.f
11. model/des/mpi_init_des_mod.f
12. model/des/mpi_pack_des_mod.f
13. model/des/mpi_unpack_des_mod.f
14. model/des/neighbour.f
15. model/des/read_res0_des.f
16. model/des/write_res0_des.f

复制后逐文件 SHA-256 校验通过 16/16，当前 code/ 与 MFIX 26.1.2 官方来源完全一致。

model/read_namelist.f 暂不复制；read_res1_des_mod.f 和 write_res1_des_mod.f 采用官方实现；cfnewvalues.f 采用官方实现；des_energy_audit.f 仅作为验证工具迁入 validation/analysis。

## 物理和数据布局决定

- 生产 JKR 只支持球形颗粒。SuperDEM/GSP 分支保持官方行为，不新增 JKR 数据、力或历史处理。
- USE_JKR=.FALSE. 是硬兼容边界：不得改变 26.1.2 原生 Hertz/LSD 的力、历史更新、墙面接触和 DMP 行为。
- PFT/PFR 是有方向的 contact-pair CC 历史。owner 从 (LL,I) 反转为 (I,LL) 时二者都反号；JKR_CONTACT_STATE 是无方向逻辑状态，不反号。
- 墙面 facet/PFT/PFR/state 是粒子附着历史，随 particle packet 迁移；颗粒–颗粒 PFT/PFR/state 属于 pair packet。
- 旧版通过 DES_USR_VAR(1,:) 搭载逐颗粒表面能。实现前必须验证其在 26.1.2 的 DMP/restart 生命周期；不能证明稳定语义时不得直接沿用。
- 邻居重建、粒子排序/所有权迁移、restart global-to-local remap 必须共享同一 global ID、CC owner 和向量方向规则。
- energy audit、历史 seed、ghost/slot/force 的 DES_USR_VAR 导出均属于验证边界，不得无条件写入生产力计算。

## Restart 格式与兼容性决定

| 来源/模式 | 决定 |
| --- | --- |
| MFIX 26.1.2 官方 DES RES 1.2 | 新 reader 可读取；缺失的 PFR/JKR/wall history 明确初始化为空，因此不能声称是活动 JKR 接触的连续续算 |
| 旧 MFIX 23.4.1 JKR DES RES 1.3 | 不直接兼容；其二维接触身份、NEIGHBOR_INDEX 和字段顺序与 26.1.2 不同。需要时另做离线转换器 |
| 新 JKR/MFIX 26.1.2 格式 | 使用独立 DES RES 2.0 标识；以官方一维 global neighbor ID/map 为基底，在版本化尾部写 PFR、pair state、wall facet/PFT/PFR/state |
| RESTART_1 | 必须恢复全部接触历史，并在首次邻居搜索和首次接触力计算前后验证 |
| RESTART_2 | 保留官方“重新初始化、不恢复接触历史”的语义，只作为敏感性对照 |

实现 DES RES 2.0 前，必须先写字段顺序表、元素类型/维度、串行与 DMP scatter 规则，以及 1.2/1.3/2.0 的显式 version gate 测试。

## 已有验证证据

| 日期 | 基准/工况 | 结果 | 结论 |
| --- | --- | --- | --- |
| 2026-08-21 | MFIX 23.4.1 旧 JKR，restart_history continuous vs 5 ns RESTART_1 | 失败 | 约 6 ns 首先出现 PFT 差异；PFR 在续算后接近零，与 23.4.1 GLOBAL_TO_LOC_COL 被注释一致 |
| 2026-08-21 | 同一工况的 RESTART_2 | 预期不一致，敏感性对照有效 | 工况确实依赖历史 |
| 2026-08-22 | MfixVer26/model vs 安装的 26.1.2 model | 一致 | 前者可作为只读官方快照 |
| 2026-08-22 | code/ 16 个官方同名覆盖文件的串行基线构建 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0；solver SHA-256 e8da3c67cf2bb1b7e99299c65b70ccdfe78a9cbf9373a792a152947d498ef3d6 |
| 2026-08-22 | MFIX 26.1.2 自定义 JKR 编译/运行 | 未执行 | 第二阶段仅编译官方内容覆盖基线，不能据此宣称 JKR、restart 或物理通过 |

现有 restart 回归只覆盖串行颗粒–颗粒、正重叠切分点；尚未覆盖负重叠颈部、墙面历史、DMP owner migration、不同 MPI 分解和双精度 restart 边界快照。VTP Float32 不能作为历史保存正确性的唯一证据。

## 最大风险

1. DMP pair history 不是“追加四个字段”即可完成：官方 26.1.2 的 pair 循环、ADD_PAIR 和 PFT pack/unpack 当前没有形成有效闭环。
2. calc_force_dem.f 和 calc_collision_wall_mod.f 已大幅重构，必须逐物理分支重写，并证明 USE_JKR=.FALSE. 零差异。
3. half-cell 4*MAX_RADIUS 修正和跨 rank 对称 pair policy 均未被确认上游解决，可能造成漏接触、重复接触或 rank-dependent history。
4. restart、neighbor rebuild、des_functions 排序/remap 和 DMP owner migration 对同一 CC 历史有多条入口，任何遗漏都会造成静默错位。
5. 旧 JKR RES 1.3 与新版官方 RES 1.2 不能按版本号大小推断兼容。
6. 旧力学文件混有生产逻辑、能量审计和验证 seed/导出，迁移时必须拆分。

## 建议迁移顺序

1. 建立项目和原生 26.1.2 baseline，固定串行 USE_JKR=.FALSE. 对照。
2. 迁移 namelist、参数检查、数据声明、分配/grow/deallocate，并立即编译。
3. 明确逐颗粒表面能字段及其 DMP/restart 生命周期，实现异质等效参数。
4. 在新版 calc_force_dem.f 上实现颗粒–颗粒 JKR，再加入 PFT/PFR/rolling/state。
5. 在新版墙面接口上实现 JKR 和 wall history，先串行验证 facet slot 连续性。
6. 扩展 neighbour.f 和 des_functions.f，验证同向/反向 owner 的历史映射。
7. 用单一协议表重写 DMP packet 的 size、pack、unpack，测试固定两进程所有权迁移。
8. 基于官方 collision restart map 实现 DES RES 2.0，完成串行后再做 DMP restart。
9. 最后加入隔离的 energy audit 和双精度边界诊断。
10. 两颗粒、墙面、restart、DMP 通过后，再运行小规模团聚体集成测试。

## 已知限制和待完成事项

- 新 Git 仓库、code/ 和 validation/{cases,analysis,results} 已建立。
- 16 个 26.1.2 生产源码已复制但内容未修改；第三阶段尚未开始。
- 串行官方内容覆盖基线已编译；SMP/DMP 和后续 JKR 修改尚未编译。
- mpifort 必须在完整激活的 mfix-26.1.2 环境中重新验证。
- 尚未决定旧 JKR 1.3 restart 是否值得开发离线转换器；默认不兼容。
- 尚未生成旧存档删除清单，也未删除、移动或压缩任何旧文件。

## 第二阶段记录

日期：2026-08-22。

- 在 /home/yuxing/mfix/MfixVer26 初始化 main 分支 Git 仓库。
- 建立扁平 code/、validation/cases、validation/analysis 和 validation/results。
- model/ 作为现有 26.1.2 官方只读快照保留，但由 .gitignore 排除。
- 创建 README.md 和 .gitignore；大型结果、solver/build、MFIX 输出、Python 缓存和临时比较文件均不进入 Git。
- code/ 的 16 个文件均直接复制自本地 model/，SHA-256 与官方来源一致。
- 本阶段已成功编译串行官方内容覆盖基线；没有运行工况或修改 Fortran，因此没有新的物理或 restart 通过结论。

## 第三阶段模块 1：参数与历史数据生命周期

日期：2026-08-22。以下状态取代上文第二阶段的“16 个文件均未修改”基线描述。

- 新增 JKR namelist：USE_JKR、颗粒/墙面表面能、颗粒/墙面临界滚动角和拉伸截止安全系数；USE_JKR 默认关闭。
- JKR 仅支持普通球形 Hertzian DEM。与 USE_COHESION、LSD、SuperDEM 或 glued-sphere 组合时启动报错。
- 新增 pair PFR、pair JKR contact state、wall PFR 和 wall JKR contact state；已覆盖 allocate、neighbor grow、particle grow、new-particle initialization、particle swap 和 deallocate。
- 新增官方 des_init_arrays.f 项目副本，仅扩展墙历史初始化；新增官方 read_namelist.f 项目副本，使其优先包含 code/desnamelist.inc。环境中的官方源码未修改。
- 本模块只建立输入、校验和状态存储，不宣称 JKR 法向力、滚动力矩、邻居重建、DMP 或 restart 已完成。

| 2026-08-22 | 模块 1 串行编译 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0、CMake 4.2.3；17 个 Fortran UDF + desnamelist.inc 被采用 |
| 2026-08-22 | USE_JKR=.FALSE. 两球 Hertz 短时 smoke | 通过 | 2 粒子完成初始化并推进；旧 VTK_PART_COL_FORCE 在 26.1.2 中已删除，测试副本移除后通过 |
| 2026-08-22 | USE_JKR=.TRUE. 关键词/分配 smoke | 通过 | 新关键词被解析，新增数组完成初始化并运行；此结果不作为 JKR 力学验证 |

下一提交将按 MFIX 26.1.2 的 Hertz 归一化契约实现逐颗粒表面能和异质接触派生参数；不会直接复制 23.4.1 的派生系数块。

## 第三阶段模块 2：逐颗粒表面能与异质参数

日期：2026-08-22。

- DES_USR_VAR(1,particle) 是生产逐颗粒表面能槽：正值覆盖 phase 默认值，非正值回退并归一化为 DES_SURFACE_ENERGY(phase)。DES_USR_VAR(2:31) 不进入生产物理。
- MFIX 26.1.2 官方路径已覆盖该粒子属性的初始化、grow、swap、DMP ghost/owner packet 和 RESTART_1 双精度读写；PFT/PFR/contact state 不借用此槽。
- 异质颗粒对使用 gamma_eff=sqrt(gamma_i*gamma_j)；phase 默认 pair/wall 的 a0、delta_c、pull-off force、刚度和 Chen damping 基系数在 Hertz 参数检查阶段初始化。
- JKR 系数采用独立的接触半径契约：运行时 KN/KT 乘 a，阻尼基系数乘 sqrt(a)；未把 23.4.1 原始 Hertz 数组布局覆盖到 26.1.2 的归一化数组。

| 2026-08-22 | 模块 2 串行编译 | 通过 | GNU Fortran 15.3.0；无未定义值或浮点异常 |
| 2026-08-22 | 同相异质表面能 0.01/0.04 | 通过 | 双精度 particle_output.csv 保持两值；几何平均等效 gamma 为 0.02 |

Restart 待办：26.1.2 read_res0_des.f 在 restart 中的 DES_USR_VAR 数量大于当前 DES_USR_VAR_SIZE 时会跳过读取而不安全地推进记录。新 reader 必须做显式容量错误检查；含逐颗粒 gamma 的 JKR restart 要求 DES_USR_VAR_SIZE>=1。

## 第三阶段模块 3：球形颗粒–颗粒 JKR、PFT 与 Chen 滚动

日期：2026-08-22。

- 在 MFIX 26.1.2 的普通球形 DEM 接触循环内语义重写 JKR 法向力；正重叠建立颈部，负重叠仅延续已有颈部，达到 `JKR_ALPHA_SAFETY*delta_c` 后清除 PFT、PFR 和 contact state。
- 接触半径由稳定分支 Newton 迭代求解；JKR 法向刚度、切向刚度和阻尼分别按 `a`、`a` 和 `sqrt(a)` 缩放。逐颗粒表面能在接触时计算异质 `gamma_eff`。
- PFT 继续沿用 26.1.2 的接触槽与切平面投影；JKR 滑移上限使用含拉脱力的弹性法向载荷。Chen PFR 按 global-ID canonical frame 累积、限幅并产生等大反向滚动力矩。
- `USE_JKR=.FALSE.`、零表面能以及 GSP/SuperDEM 路径均不进入 JKR 分支；邻居重建、DMP 和 restart 中的 PFR/state 搬运仍由后续模块完成。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | 模块 3 串行编译 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0；`git diff --check` 通过 |
| 2026-08-22 | `USE_JKR=.FALSE.` 两球 Hertz 对照 | 通过 | 修改前后 `particle_output.csv` SHA-256 均为 `7cdc2bf6a09e68480dd2ebd52183b0ca50ae02c4326738ebb57155f5a9d4f8b6` |
| 2026-08-22 | `pair_normal_stick` | 通过 smoke | 139758 个 DEM 子步完成，无 JKR Newton 警告；末态保持对称运动。尚需双精度接触状态诊断确认负重叠颈部全过程 |
| 2026-08-22 | `pair_normal_rebound` | 通过 smoke | 两粒子末态速度为 `-0.11228027/+0.11228027 m/s`，形成拉脱分支 |
| 2026-08-22 | `pair_oblique_stick` | 通过 smoke | 末态存在切向速度与约 `2.5367e5 rad/s` 对称自旋，证明切向/滚动路径被激活；PFT/PFR 数值仍需 validation-only 双精度快照 |

本模块的 smoke test 只证明执行路径和宏观响应；它不替代后续力–重叠、接触半径、PFT/PFR、contact state 和解析拉脱力的定量回归。

## 第三阶段模块 4：球形颗粒–墙面 JKR 与墙历史

日期：2026-08-22。

- 基于 26.1.2 的 STL facet 搜索、法向构造和力/力矩累加接口实现球形墙面 JKR；未进入 SuperDEM 或 glued-sphere 路径。
- JKR 候选范围扩展到 `delta_c`，实际拉伸截止使用 `JKR_ALPHA_SAFETY*delta_c`；正重叠建立墙颈部，负重叠只延续已有 facet-slot state。
- 墙面 PFT 保留官方切向历史算法；新增 Chen PFR、临界滚动位移和 `NORMAL x PFR` 滚动力矩。
- facet slot 的新建、陈旧 slot 复用、显式移除、无邻近 facet 和动态扩容均同步初始化或清除 facet ID、PFT、PFR 和 wall JKR state，防止历史串槽。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | 模块 4 串行编译 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0；`git diff --check` 通过 |
| 2026-08-22 | 官方 26.1.2 vs 自定义 solver 墙面 Hertz | 逐字节一致 | 两份 `particle_output.csv` SHA-256 均为 `b560bb90ef4221b2c0799a64a59b0e0a4f83a33c295507c774d8391158284bc9` |
| 2026-08-22 | `wall_normal_capture` | 通过 smoke | 单粒子完成 1 us，末态中心–墙间隙约 3.70 nm，处于负重叠颈部范围；无 Newton 警告 |
| 2026-08-22 | `wall_oblique_stick` | 通过 smoke | 完成 0.6 us，末态切向速度约 0.7793 m/s、自旋约 -7.835e5 rad/s；切向/滚动路径激活 |

墙面 restart、DMP owner migration 和双精度 wall facet/PFT/PFR/state 快照尚未完成，因此本模块只提交串行生产物理与内存生命周期。

## 第三阶段模块 5：邻居重建与粒子排序历史

日期：2026-08-22。

- 在官方 26.1.2 neighbor rebuild 的 PFT 快照/恢复上同步追加 PFR 和 JKR contact state；同向 pair 原值恢复，新 pair 清零/false。
- JKR 模式增加反向 owner fallback：PFT/PFR 取反，标量 contact state 不取反；该 fallback 在 `USE_JKR=.FALSE.` 时不启用。
- 普通球形 JKR 粒子排序使用独立的双精度历史例程，同时搬运 PFT、PFR 和 state；官方排序例程与 SuperDEM 路径保持不变，以维持关闭 JKR 时的原生行为。
- 当前反向 fallback 仍依赖重排后的 local particle index；DMP owner migration 必须在 packet 层使用 global ID pair 规则补全。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | 模块 5 串行编译 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0 |
| 2026-08-22 | JKR `pair_oblique_stick` 改动前后 | 逐字节一致 | `particle_output.csv` SHA-256 均为 `b5681e59fe49c562a1bcd82cc1baa51e4045d8252f0c3c9d32aa991ae6fdbce8` |
| 2026-08-22 | 官方墙 Hertz vs 模块 5 后 solver | 逐字节一致 | SHA-256 仍为 `b560bb90ef4221b2c0799a64a59b0e0a4f83a33c295507c774d8391158284bc9` |

尚需 validation-only 双精度快照直接证明邻居重建前后 PFT/PFR/state 连续；本模块不宣称 DMP 或 restart 已通过。

## 第三阶段模块 6：DMP 所有权迁移

日期：2026-08-22。

- 基于 26.1.2 的 particle-owner packet 追加墙面 facet ID、PFT、PFR 和 JKR state；基于现有 11-field pair packet 发送 global owner/neighbor ID、网格位置、PFT、PFR 和 pair state。
- pair history 不在 unpack 时猜测 CC，而是按 global-ID pair 暂存；官方邻居重建产生有效 CC 后再恢复。同向向量保持，反向 owner 映射时 PFT/PFR 取反，state 不取反。
- `USE_JKR=.FALSE.` 时不发送新增 pair/wall 历史，保持官方 26.1.2 owner migration 路径。
- DMP 粒子包容量在初始化时由墙面槽数固定。为避免某一 rank 动态扩容后产生不一致包布局，JKR DMP 在 8 个墙面 facet slot 全部耗尽时显式停止；串行仍保留官方动态扩容。当前适用范围确认不会超过默认 8 槽，因此不要求扩展动态容量。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | DMP solver 编译 | 通过 | MFIX 26.1.2、GNU Fortran 15.3.0、MPI 3.1；`ENABLE_MPI=ON` |
| 2026-08-22 | `pair_new_contact_dmp`，1 rank vs 2 ranks | 通过 smoke | 两者均完成 0.65 us；4 粒子最终 X/Y/Z/U/V/W 双精度 CSV 逐项差值为 0 |
| 2026-08-22 | `wall_cross_dmp`，2 ranks | 通过 smoke | 单粒子完成 0.15 us；日志无 error、warning、abort、NaN 或 overflow |

临时运行目录为 `/tmp/jkr26_dmp_pair_np1`、`/tmp/jkr26_dmp_pair_np2` 和 `/tmp/jkr26_dmp_wall_np2`，不进入 Git。

上述结果验证 DMP 可运行及宏观串并行一致性，但尚未在所有权交换的四个指定断点直接导出 PFT/PFR/JKR state。因此“跨分区历史逐字段保持”仍列为分层验证待办，不以本 smoke test 代替。

## 第三阶段模块 7：RESTART_1 历史格式

日期：2026-08-22。

- `USE_JKR=.FALSE.` 继续写 MFIX 26.1.2 官方 DES restart 1.2；`USE_JKR=.TRUE.` 写项目格式 2.0。2.0 保留完整官方 1.2 前缀，再追加 pair owner/neighbor global ID、PFR、JKR state，以及墙面 facet ID/PFT/PFR/state。
- reader 先采用官方 `READ_PAR_COL` 完成 collision ownership 分配；由于官方 reader 不恢复旧 `NEIGHBOR_INDEX`，2.0 再按 global-ID pair 放入与 DMP 共用的 pending queue，首次官方邻居重建后恢复到有效 CC。墙历史通过 PARRAY 继续按 particle global ID 映射。
- 旧 23.4.1 JKR 1.3 与新布局不兼容，reader 显式拒绝 1.3；JKR RESTART_1 不允许读取缺少历史的 1.2。RESTART_2 保留官方碰撞历史前 early return，可作为历史敏感性对照。
- reader 对 `DES_USR_VAR_SIZE` 不足改为显式错误，避免官方旧逻辑跳过记录后继续错位读取。新 2.0 不承诺读取旧 JKR RES；兼容决定以明确失败代替静默误读。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | 串行与 DMP 编译 | 通过 | 最新产物分别为 `code/build/mfixsolver` 与 `code/build_dmp/mfixsolver_dmp` |
| 2026-08-22 | 等径两球 continuous vs 5 ns + RESTART_1 | 逐字节一致 | 双精度 CSV SHA-256 均为 `5b2823296d83d0471bbbfdc2827a8e8555856b24ff10e985a77dba85e53e2831` |
| 2026-08-22 | 异径两球非零 PFR restart | 逐字节一致 | 5 ns 读回 PFT=`8.17869e-10 m`、PFR=`7.88762e-12 m`、state=TRUE；末态 SHA-256 均为 `8cf9c6fac93f7bce696d4d0ab56c569bd38f8d065bf3b8eda959f47fea32f309` |
| 2026-08-22 | 墙面斜碰 continuous vs 0.2 us + RESTART_1 | 逐字节一致 | 双精度 CSV SHA-256 均为 `b3ad295807141728ce7580fda695aba6e261d1706edaa57a63ab29b51ba0f402` |
| 2026-08-22 | RESTART_2 敏感性 | 通过对照 | 末态 SHA-256 `83a18505d2afc06b95743faaf48afb19ee5964f2e871cda89a8f01b61d1a4714`，与 continuous 不同 |
| 2026-08-22 | 最终 `USE_JKR=.FALSE.` 墙 Hertz | 逐字节一致 | SHA-256 保持 `b560bb90ef4221b2c0799a64a59b0e0a4f83a33c295507c774d8391158284bc9` |
| 2026-08-22 | 固定 2-rank continuous vs 2-rank RESTART_1 | 暂未通过 | 两者均正常完成，但末态 V 分量最大差约 `1.18e-4 m/s`；需用双精度断点快照区分官方 DMP restart 重排差异与 JKR queue 时序 |
| 2026-08-22 | 串行写入、2-rank 读取 | 暂未通过 | 可完成但与串行 continuous 有可见偏差；不同分解兼容不列入当前通过项 |

串行 RESTART_1 已覆盖 PFT、非零 PFR、pair neck state 和墙面历史；DMP restart 仍是验收阻塞项，不因运行完成而标记通过。

## 第三阶段模块 8：粗网格邻居搜索修正与 DMP restart 基线

日期：2026-08-22。

- MFIX 26.1.2 官方 `desgrid_mod.f` 的 half-cell coarse 判据仍使用 `2*MAX_RADIUS*FACTOR_RLM`。旧 JKR 修改的意图是比较半个网格宽度与最大搜索半径，因此在 26.1.2 邻居框架上语义重写为系数 4；没有复制旧文件或改变新版其余 DES grid/DMP 过滤逻辑。
- 粗网格新接触工况在串行与固定两进程 DMP 下均形成一个激活的 JKR contact state，最终 4 个颗粒的 X/Y/Z/U/V/W 双精度 CSV 逐项差值为 0；该工况为法向接触，PFT/PFR 预期为零。
- DMP RESTART_1 首次邻居重建后的 validation-only 快照直接观察到拥有接触对的 rank。等径 2-rank 写/读工况恢复 PFT=`8.2190147440577080e-10 m` 与 state=TRUE；异径串行写档、2-rank 读取工况恢复 PFT=`8.1786908176533694e-10 m`、非零 PFR=`7.8876185564904887e-12 m` 与 state=TRUE。临时输出均已从生产源码撤除。
- 固定两进程 `USE_JKR=.FALSE.` 官方 1.2 restart 基线本身也不是 continuous 逐字节一致：末态最大差 X=`7.9490e-11 m`、Y=`1.3461e-13 m`、U=`2.5411e-4 m/s`、V=`2.7856e-6 m/s`，日志无错误。因此模块 7 中 JKR DMP restart 的非零末态差异不能全部归因于新增历史，但仍未达到本项目 DMP restart 验收标准。

| 日期 | 工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | `pair_new_contact_coarse_dmp`，1 rank vs 2 ranks | 通过 | 均退出 0；4 粒子末态 X/Y/Z/U/V/W 最大差均为 0；形成 1 个 JKR 接触 |
| 2026-08-22 | JKR RESTART_1 首邻居重建，2-rank reader | 逐字段恢复通过 | 非零 PFT、非零 PFR 与 state 均按 global-ID pending queue 恢复；末态 continuous 一致性仍受官方 DMP restart 基线差异阻塞 |
| 2026-08-22 | 2-rank Hertz 官方 1.2 continuous vs RESTART_1 | 基线非零差异 | 无运行错误；用于隔离 MFIX 26.1.2 官方 DMP restart 重排误差，不作为 JKR 通过结论 |

跨 rank 单侧 pair 过滤策略尚未被独立证明在所有 owner 迁移方向都正确；后续仍需用非零 PFT/PFR 的跨分区接触逐断点检查，不能由本模块的法向 coarse case 替代。

## 第四阶段：固化回归与团聚体集成 smoke

日期：2026-08-22。

- 新增 `validation/cases/restart_pair_pfr/`：异径两球、非零 PFT/PFR 和 JKR neck state 的 continuous、5 ns + RESTART_1、5 ns + RESTART_2 runner。case 仅保留生产 `DES_USR_VAR(1)` 表面能槽，2:31 能量/调试槽未混入生产源码。
- 新增 `validation/cases/restart_wall_oblique/`：斜碰墙面的 continuous、0.2 us + RESTART_1 runner，覆盖 wall facet/PFT/PFR/state 的集成恢复。
- 两个 runner 只清理各自在 `validation/results/` 下的固定目录，接受相对或绝对 solver 路径；大型输出继续由 `.gitignore` 排除。
- 复用旧 `initialAgg.mfx` 与 1000 粒子输入，在 MFIX 26.1.2 串行 solver 下推进一个有效宏步至 10 us。旧输入 `Diameter: F` 的零默认值仅在 `/tmp` 副本修正为原材料直径 2 um，未修改旧存档。

| 日期 | 固化/集成工况 | 结果 | 证据/边界 |
| --- | --- | --- | --- |
| 2026-08-22 | `restart_pair_pfr` runner | 通过 | continuous 与 RESTART_1 SHA-256 均为 `b80cf63a7b052c198fc5c632ab75bc2b50df5c5e7033aabefaabe5bdf454806d`；RESTART_2 为 `3142dea040e6d3cc704836250f27d95fc3a60d474790c3f51250bc8e3946e980` |
| 2026-08-22 | `restart_wall_oblique` runner | 通过 | continuous 与 RESTART_1 SHA-256 均为 `b3ad295807141728ce7580fda695aba6e261d1706edaa57a63ab29b51ba0f402` |
| 2026-08-22 | 1000 粒子 initialAgg，一个宏步 | 通过 smoke | 达到 10 us，1000/1000 保留；Rg=`9.6850117176279204e-6 m`，平均速率=`1.4091528096525046e-2 m/s`；日志无 error/warning/abort/NaN/overflow |

团聚体结果只证明输入兼容、数组/接触路径可运行及短时宏观量有限；尚未计算平均配位数、重叠率、去除向心力后的长期稳定性或 restart 宏观连续性，不能替代两颗粒历史回归。

## 当前验收分类

| 分类 | 当前决定 |
| --- | --- |
| 已迁移 | 球形颗粒–颗粒/墙面 Chen/JKR，负重叠颈部，PFT/PFR/滚动力矩，逐颗粒表面能，邻居重建，DMP owner packet，DES RES 2.0 RESTART_1 pair/wall history，以及固定 2-rank/2×1×1 DMP restart 验收 |
| 采用新版官方实现 | MFIX 26.1.2 `READ_PAR_COL`/global-ID collision map、RESTART_2 early return、粒子所有权与 ghost 框架、SuperDEM/GSP 原生代码路径；在其上仅追加 JKR 历史 |
| 暂缓/范围外 | 通用不同 MPI 分解遍历、所有方向的跨 rank 反向 owner 断点、超过默认 8 个 wall facet slots、完整能量审计、团聚体长期/配位/重叠/restart 指标 |
| 废弃 | 旧 23.4.1 JKR RES 1.3 直接兼容、旧式 `NEIGHBOR_INDEX` restart、旧 `read_res1` 缺陷修补、`cfnewvalues.f` 能量审计调用、生产代码中的 DES_USR_VAR(2:31) seed/诊断、旧 `read_namelist.f` 空白差异 |

`USE_JKR=.FALSE.` 的串行 Hertz/墙面基准保持 MFIX 26.1.2 原生逐字节一致。DMP restart 不要求绕过 MFIX 26.1.2 自身的并行重排差异，而采用双重判据：JKR 边界历史必须精确相等，最终轨迹必须满足固定容差且不得差于同工况 JKR-off Hertz 基线加绝对容差。

## 第五阶段：固定两进程 DMP restart 验收

日期：2026-08-24。

- `desgrid_mod.f` 的 `4*MAX_RADIUS*FACTOR_RLM` half-cell 判据现仅在 `USE_JKR=.TRUE.` 时启用；JKR-off 路径恢复 MFIX 26.1.2 官方 `2*MAX_RADIUS`，避免把 JKR 搜索修正扩散到原生 Hertz。
- 在 pair 首次力计算前、wall restart 读回和 wall 首次力计算前补齐通用双精度诊断事件。wall 读回诊断只遍历与力循环一致的 `IS_NORMAL` 粒子，避免把 `MAX_PIP` 预分配空槽的未定义值误认作接触。诊断由 `MFIX_JKR_DIAGNOSTICS` 显式启用，初始化和文件写入受 OpenMP critical 保护；默认运行不创建诊断文件。
- `dmp_restart` 对颗粒对和斜碰墙面分别运行 JKR 与 JKR-off Hertz 的 continuous、stage1、RESTART_1、RESTART_2，共 16 个两进程求解分支。VTP 仅用于按 global particle ID 对齐最终轨迹；PFT、PFR 和 JKR state 使用双精度内部事件精确比较。
- wall history 保持 MFIX/当前迁移的默认固定 8 槽。当前适用范围明确假定单颗粒不会同时接触超过 8 个壁面，不为未预期工况引入动态 packet 布局。

| 工况 | 边界历史 | JKR 最终态 | 相对 Hertz 基线 | 结果 |
| --- | --- | --- | --- | --- |
| `restart_pair_pfr`，2 ranks/2×1×1 | 1 个 global-ID pair 从 `unpack_queue`、pending restore、neighbor restore 到首次 force 的 PFT/PFR/state 精确相等 | 固定 DMP VTP 容差内；最大位置误差 `9.09e-13 m`、速度 `2.09e-5 m/s`、角速度 `28.60 rad/s`、碰撞力 `1.99e-12 N` | 全字段小于 Hertz 基线加绝对容差 | PASS |
| `restart_wall_oblique`，2 ranks/2×1×1 | 1 个 `(particle_gid,facet_id)` 从 restart read 到首次 wall force 的 PFT/PFR/state 精确相等 | 四个比较字段末态误差均为 0 | 不差于零误差 Hertz 基线 | PASS |

串行和 DMP solver 均重新编译通过；统一 runner 中 99 项单元测试全部通过。当前 solver SHA-256 分别为 `09feaa7786e257b932b2e1a022bf449c47be089ee55c28e12cd21a6e45502e70` 和 `6e1a05241cee6c301cb405f93f72bd99827abe1277a66d5add060b9df1576fb7`。诊断钩子属于为可重复验收保留的默认关闭、非物理改动；按当前决定暂不迁移。

## 仓库结构迁移

日期：2026-08-24。

- 生产仓库继续只跟踪 `code/`、`docs/`、根 `README.md` 和 `.gitignore`；`model/`、`code/build*` 与嵌套的 `validation/` 均保留在工作区但不进入生产 Git。
- 原 `MfixVer26TestCases` 仓库连同 `.git` 整体移动到 `MfixVer26/validation`，再以普通 rename commit 将原内层 `validation/*` 扁平化到 Validation 仓库根。旧 commit ID 和历史保持不变。
- Validation 仓库直接包含 `cases/`、`jkr_tools/`、`tests/`、`results/` 和运行脚本；`results/` 中仅 `.gitkeep` 被跟踪，计算结果继续保留但忽略。
- 迁移前两个仓库的完整 bundle 保存于 `/home/yuxing/mfix/git-backups/2026-08-24-validation-layout/`，并已通过 `git bundle verify`。

## 第六阶段范围确认

日期：2026-08-24。

- 已完成的异质表面能、固定 2-rank DMP 和 restart 边界计划归档为
  `docs/archive/NEXT_STAGE_PROMPT_2026-08-23.md`；当前计划继续写在
  `docs/NEXT_STAGE_PROMPT.md`，避免历史完成项与新待办混用。
- 默认关闭的 `MFIX_JKR_DIAGNOSTICS` 暂不迁移。已有法向 capture/rebound
  与独立 JKR 参考覆盖负重叠颈部及临界断裂，不重复建立 tensile-neck case。
- 已确定的实施顺序为定量斜碰/滑动/滚动、wall face/edge/corner/seam、以及
  2×2×2 双尺度颗粒群 serial–DMP 对照；三项现均已完成。最后一项只比较
  同模型的串并行结果：
  JKR 对 JKR、JKR-off 原生 MFIX 对 JKR-off，不要求两种物理模型彼此一致。

## 第六阶段任务 1：定量斜碰、滑动与滚动

日期：2026-08-24。

- Validation 仓库新增 `tangential` suite 和两个 canonical case，分别覆盖
  不等径颗粒对和颗粒–墙面的斜向运动、stick/slip、PFT/PFR、滚动饱和及
  非零角速度响应；生产 Fortran、默认参数和默认关闭诊断均未修改。
- 独立分析复现 `KT=8G*a`、Chen damping、JKR Coulomb 上限、滚动临界
  位移与滚动刚度。不等径颗粒对的 damping 明确按 MFIX 相级参考质量
  `REF_MASS/2` 计算；若误用实际粒子约化质量，会产生约 8% 的伪差异。
- 颗粒对专用 case 仅把输入 `jkr_thetar_crit` 设为 `0.001`，以便在短时
  受控工况中触发滚动上限；这不是生产默认值变更。

| 工况 | stick/slip events | 最大切向力残差 | Coulomb/滚动最大比 | 结果 |
| --- | ---: | ---: | ---: | --- |
| `pair_oblique_sliding_rolling` | 3489/613 | `7.04e-21 N` | `1.0000000000000016 / 1.0000000000000004` | PASS |
| `wall_oblique_sliding_rolling` | 81608/2247 | `1.53e-19 N` | `1.0000000000000018 / 1.0` | PASS |

`python3 -m unittest discover -s tests -v` 共运行 106 项，其中 101 项通过，
5 项因缺少可选 `matplotlib`/`fparser` 依赖跳过。下一项转入 wall
face/edge/corner/triangulated seam。

## 第六阶段任务 2：真实 STL 壁面几何

日期：2026-08-24。

- Validation 生成器新增可审查的 ASCII `geometry_0001.stl` 资产及
  `wall_geometry` suite。single face、正交 edge、三面 corner 和两块共面
  三角形 seam 均使用 `cartesian_grid=.True.`、`use_stl=.True.`，不再沿用
  MFIX23 中实际为 box boundary 的同名 smoke case。
- 初次 seam 运行在零重叠处同时报告 facet 1/2，两者标量 JKR 力完全相同，
  VTP 合力约为独立单面理论值的两倍；face、edge、corner 分别正确报告
  1、2、3 个物理贡献。这证明旧测试无法发现真实 STL 共面双计。
- `calc_collision_wall_mod.f` 现仅在有效 JKR 分支中，对法向一致、最近点
  相同且已有活动 JKR state 的共面 facet 保留第一项。正交法向不会合并，
  `USE_JKR=.FALSE.` 不进入该逻辑，因此没有改变 MFIX26 原生 Hertz 路径。

| 工况 | 活动 facet ID | 双精度标量力最大误差 | VTP 合力相对误差 | 结果 |
| --- | --- | ---: | ---: | --- |
| `wall_stl_single_face` | `[1]` | `3.02e-21 N` | `9.35e-5` | PASS |
| `wall_stl_edge` | `[1,2]` | `3.02e-21 N` | `1.87e-4` | PASS |
| `wall_stl_corner` | `[1,2,3]` | `3.02e-21 N` | `2.81e-4` | PASS |
| `wall_stl_facet_seam` | `[1]` | `3.02e-21 N` | `9.35e-5` | PASS |

串行和 DMP solver 均编译通过；110 项单元测试中 105 项通过、5 项因可选
依赖缺失跳过；7 个官方 MFIX26 与迁移 solver 的 JKR-off Hertz 对照全部
PASS。新 solver SHA-256 分别为
`46c16ba38f63af4902e14f4c1e70403711a9d5f3ac7a0cb977c95ff2ac1fcbb6`
和 `6578c5625953fc573879afd8ddecdf7fa12e89b172b055c067a0cd4418c74de1`。
随后转入固定 2×2×2、8-rank 双尺度颗粒群 serial–DMP 对照，结果见下节。

## 第六阶段任务 3：固定 8-rank 双尺度颗粒群

日期：2026-08-24。

- Validation 新增 `cloud_8rank` suite 和四个 canonical case，固定
  8 ranks/2×2×2，不扫描其他 rank 或分解。每个 rank 初始一个粒子，四组
  pair 跨 x 分区相向运动；serial 与 DMP 除 `nodesi/j/k` 外物理输入字节一致。
- `R=1 μm`、`FACTOR_RLM=1.2` 时 DES search range 为 `2.4 μm`，两倍为
  `4.8 μm`。低尺度 `W=3.333 μm` 位于 MFIX 官方系数 2 与 JKR 系数 4
  判据之间；高尺度 `W=5.0 μm` 位于两者之上。
- 低尺度采用确定性 half-cell 反例排列并交替 pair 的 global-ID 顺序。
  JKR 的系数 4 路径保留四组接触且 serial/DMP 全字段相同；JKR-off 原始
  系数 2 路径稳定复现 rank-dependent 接触集合差异。这一分支被显式标记为
  负对照，PASS 表示预期缺陷出现，不表示原始路径一致。

| 尺度/模型 | serial–DMP 结果 | 接触集合 | 最大绝对误差 x/v/omega/force | 验收 |
| --- | --- | --- | --- | --- |
| `3.333 μm` / JKR | 一致 | 四对相同 | `0 / 0 / 0 / 0` | PASS |
| `3.333 μm` / JKR-off | 预期不一致 | `2.4000e-7 s`：serial 2 对、DMP 4 对 | `4.36e-7 / 1.0301 / 0 / 2.10e-7` | PASS（负对照） |
| `5.0 μm` / JKR | 一致 | 四对相同 | `0 / 0 / 0 / 0` | PASS |
| `5.0 μm` / JKR-off | 一致 | 四对相同 | `0 / 0 / 0 / 0` | PASS |

JKR 两档和高尺度 JKR-off 均比较 101 帧、每字段 2424 个值，按 global ID
对齐后的首次物理差异为空；初始 owner 覆盖 rank 0–7。四个工况日志无
error、warning、abort、NaN 或 overflow。该专项为纯法向搜索判据工况，
PFT/PFR 预期为零；非零切向历史的跨 owner 连续性由既有
`pair_active_contact_migration` 和 DMP restart 双精度证据承担。

## 第六阶段任务 4：最终综合回归

日期：2026-08-24。

- 从两个干净仓库执行 `./run_validation.sh all`，统一结果为 42/42 PASS；
  精确来源 commit 由最终快照的 `RUN_INFO.txt` 记录。
- 117 项 unit tests 中 113 PASS，4 项因缺少可选 `fparser` 跳过。数值回归
  覆盖官方 Hertz 保护、serial JKR、定量切向、true-STL、2-rank DMP、
  8-rank 双尺度、serial restart 和 DMP restart。
- 77 份 solver 日志无 error、abort、NaN 或 overflow。保留并明确记录
  4 条官方 wrapper 的 namelist description 截断提示，以及 single-face/
  seam STL 各 1 条 `UNABLE TO PROPAGATE F_AT ARRAY` warning；对应官方轨迹
  与 STL facet/力/state 验收均通过，不将提示隐藏为“零 warning”。
- 最终 solver SHA-256 为 serial
  `46c16ba38f63af4902e14f4c1e70403711a9d5f3ac7a0cb977c95ff2ac1fcbb6`
  和 DMP
  `6578c5625953fc573879afd8ddecdf7fa12e89b172b055c067a0cd4418c74de1`。
