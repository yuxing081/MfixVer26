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

所有命令均在 WSL Ubuntu-22.04 中执行。先激活准确环境：

    conda activate /home/yuxing/miniforge3/envs/mfix-26.1.2

已核验 MFIX 26.1.2、GNU Fortran 15.3.0、CMake 4.2.3 和 MPI 3.1。串行与 DMP solver 的正式产物分别位于 `code/build/mfixsolver` 和 `code/build_dmp/mfixsolver_dmp`；不要使用 `code/` 顶层可能残留的旧 solver 副本。

从干净构建目录配置并编译串行 solver：

    cd /home/yuxing/mfix/MfixVer26/code
    cmake -S "$CONDA_PREFIX/share/mfix/src" -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j 8

配置并编译 DMP solver：

    cd /home/yuxing/mfix/MfixVer26/code
    cmake -S "$CONDA_PREFIX/share/mfix/src" -B build_dmp -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON
    cmake --build build_dmp -j 8

build/、solver、模块文件和运行输出均被 Git 忽略。不要在 conda 官方源码目录内构建或编辑。

## 最小运行与测试

运行固化的异径两球 PFT/PFR/JKR state restart 回归：

    cd /home/yuxing/mfix/MfixVer26
    validation/cases/restart_pair_pfr/run_restart_pair_pfr.sh code/build/mfixsolve

脚本在 `validation/results/restart_pair_pfr/` 生成三条分支并比较双精度 `particle_output.csv`：

墙面 PFT/PFR/state 回归使用：


- continuous；
- 5 ns + RESTART_1；
- 5 ns + RESTART_2 敏感性对照。

两进程运行的基本命令为：

    mpirun -np 2 /home/yuxing/mfix/MfixVer26/code/build_dmp/mfixsolver_dmp -f mfix.dat

截至 2026-08-22，串行 pair/wall JKR、异质表面能、PFT/PFR/state RESTART_1、DMP owner migration 和 coarse-grid 新接触已通过。DMP 首次邻居重建已直接恢复非零 PFT/PFR/state；但 continuous 与 DMP RESTART_1 的末态仍非逐字节一致，且 `USE_JKR=.FALSE.` 的 MFIX 26.1.2 官方 1.2 DMP restart 也存在基线差异。该项在 MIGRATION.md 中保持“暂缓/未验收”，不能由编译或运行完成替代。

验证顺序固定为：

1. 原生 26.1.2 与 USE_JKR=.FALSE. 自定义 solver 对照。
2. 串行两颗粒和墙面单接触。
3. PFT/PFR/JKR state 的 RESTART_1 回归，RESTART_2 仅作敏感性对照。
4. 固定两进程 DMP owner migration 和 restart。
5. 小规模团聚体集成测试。

## 生产代码与验证代码边界

生产 code/ 保留 JKR 力学、PFT/PFR、JKR contact state、逐颗粒表面能、邻居/DMP/restart 生命周期所必需的逻辑。新增历史仿照 MFIX 官方 PFT 的保存、全局 ID 映射和恢复流程。

能量统计、DES_USR_VAR(2:31)、人工 PFR seed、ghost/slot/force 导出和双精度边界快照只允许存在于 validation 的具体 case、analysis 工具或显式验证 patch 中，不得无条件进入生产力计算。
