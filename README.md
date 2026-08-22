# JKR for MFIX 26.1.2

## 项目目标

本项目把基于 MFIX 23.4.1 开发的球形颗粒 Chen/JKR 接触模型迁移到 MFIX 26.1.2。迁移以新版官方源码和官方 restart 邻居映射为基底，不覆盖环境源码，也不移植旧版已知 restart 缺陷。

当前范围只包括球形颗粒。SuperDEM/GSP 异形颗粒不做 JKR 适配；其官方代码路径保持原样。

## 基准与来源

- MFIX 基准：26.1.2。
- 官方安装源码：/home/yuxing/miniforge3/envs/mfix-26.1.2/share/mfix/src。
- 本地只读快照：/home/yuxing/mfix/MfixVer26/model。
- 旧 JKR 仓库：/home/yuxing/mfix/P3/JKR/versions/code。
- 旧 JKR 提交：ea55027498a4b939a2982cbb77d137201b034ffc。
- 旧工作树状态：dirty；未提交内容仅为 restart 验证 seed 和诊断字段，不能进入生产物理。
- 迁移开始日期：2026-08-22。

每个 code/ 文件均从 MFIX 26.1.2 官方版本复制，再按迁移矩阵做最小增量修改。MFIX CMake 按 basename 用当前目录顶层的同名 Fortran 文件替换官方文件，因此 code/ 保持扁平布局。

## 目录

- code/：相对于 MFIX 26.1.2 修改的生产 Fortran 源码；不放验证 seed 或 case 专用诊断。
- validation/cases/：两颗粒、墙面、restart、DMP 和最终团聚体工况；case 专用 DES_USR_VAR 诊断放在对应 case 内。
- validation/analysis/：比较器、双精度边界快照读取和可选能量审计。
- validation/results/：可重新生成的大型输出，默认不进入 Git。
- MIGRATION.md：逐文件迁移决定、restart 格式、测试结果、限制和待办。
- model/：当前已有的 MFIX 26.1.2 官方只读快照，不纳入本仓库版本控制。

## 编译环境

先进入 WSL 并激活准确环境：

    conda activate /home/yuxing/miniforge3/envs/mfix-26.1.2

已核验 GNU Fortran 15.3.0、CMake 4.2.3 和 build_mfixsolver。DMP 构建前还必须在完整激活环境中验证 mpifort。

编译官方 26.1.2 串行基准：

    cd /home/yuxing/mfix/MfixVer26
    mkdir -p build/native-serial
    cd build/native-serial
    build_mfixsolver --batch -j 8

编译 code/ 中同名源码覆盖后的串行 solver：

    cd /home/yuxing/mfix/MfixVer26/code
    build_mfixsolver --batch -j 8

编译 DMP solver：

    cd /home/yuxing/mfix/MfixVer26/code
    build_mfixsolver --batch --dmp -j 8

2026-08-22 已验证 code/ 串行基线构建成功，16 个本地覆盖文件均被识别；生成 solver 的 SHA-256 为 e8da3c67cf2bb1b7e99299c65b70ccdfe78a9cbf9373a792a152947d498ef3d6。此结果只证明源码布局和工具链可用。

build/、solver、模块文件和运行输出均被 Git 忽略。不要在 conda 官方源码目录内构建或编辑。

## 最小运行与测试

第二阶段只建立可追踪的 26.1.2 源码基线，尚未提供已通过的 JKR case。第三阶段从 USE_JKR=.FALSE. 原生对照和串行两颗粒 smoke case 开始；在相应 case 和测试脚本实际加入后，本节记录精确命令。当前不得把“能够编译”报告为物理或 restart 通过。

验证顺序固定为：

1. 原生 26.1.2 与 USE_JKR=.FALSE. 自定义 solver 对照。
2. 串行两颗粒和墙面单接触。
3. PFT/PFR/JKR state 的 RESTART_1 回归，RESTART_2 仅作敏感性对照。
4. 固定两进程 DMP owner migration 和 restart。
5. 小规模团聚体集成测试。

## 生产代码与验证代码边界

生产 code/ 保留 JKR 力学、PFT/PFR、JKR contact state、逐颗粒表面能、邻居/DMP/restart 生命周期所必需的逻辑。新增历史仿照 MFIX 官方 PFT 的保存、全局 ID 映射和恢复流程。

能量统计、DES_USR_VAR(2:31)、人工 PFR seed、ghost/slot/force 导出和双精度边界快照只允许存在于 validation 的具体 case、analysis 工具或显式验证 patch 中，不得无条件进入生产力计算。
