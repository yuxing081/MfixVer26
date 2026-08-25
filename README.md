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
- validation/：嵌套的独立 Validation Git 仓库，拥有 cases、生成器、理论参考、测试、结果契约和验证文档；外层生产仓库整体忽略该目录。
- 本仓库不跟踪 canonical 输入、验证脚本或运行结果，仅保留生产源码及其迁移说明。
- docs/JKR_IMPLEMENTATION_SUMMARY.md：JKR 公式、变量、生产代码行号、MFIX 生命周期扩展和验证证据的精炼汇总。
- docs/MIGRATION.md：生产源码的逐文件迁移决定与实现记录。
- model/：当前已有的 MFIX 26.1.2 官方只读快照，不纳入本仓库版本控制。

旧版 `versions/test` 只用于审查测试过程和公式来源，不复制 MFIX 23 结果作为理论真值。

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

验证只有一个 shell 入口；Python 运行器负责生成输入、创建独立结果目录、记录 Git/solver/编译器/参数元数据，并在任一必需测试失败时返回非零退出码：

    cd /home/yuxing/mfix/MfixVer26/validation
    ./run_validation.sh unit
    ./run_validation.sh hertz
    ./run_validation.sh serial
    ./run_validation.sh dmp
    ./run_validation.sh restart
    ./run_validation.sh dmp_restart
    ./run_validation.sh all

当前已接入 `unit`、官方 stock 对照的串行 `hertz`、基础 JKR normal `serial`、固定两进程 `dmp`、串行 `restart` 和两进程 `dmp_restart`；未运行项目不会标记为通过。

restart 套件对两颗粒和斜碰墙面分别运行：

1. continuous；
2. split stage + RESTART_1；
3. split stage + RESTART_2 敏感性负对照。

串行 RESTART_1 必须与 continuous 的双精度 `particle_output.csv` 逐字节一致。DMP RESTART_1 还要求边界 PFT/PFR/JKR state 双精度逐字段相等、global-ID 最终轨迹满足固定容差，并且误差不超过同工况 JKR-off Hertz 基线加绝对容差。RESTART_2 是否不同单独记录，不作为正向通过条件。每次 runner 调用替换 `validation/results/` 快照。

完整测试顺序保持为原生 Hertz 保护、串行墙面/颗粒对、速度扫描、壁面几何、2-rank DMP、RESTART_1 和小颗粒群集成测试。正式结果写入 `docs/VALIDATION.md`，活动进度写入 `docs/TEST_PROGRESS.md`。

## 生产代码与验证代码边界

生产 code/ 保留 JKR 力学、PFT/PFR、JKR contact state、逐颗粒表面能、邻居/DMP/restart 生命周期所必需的逻辑。新增历史仿照 MFIX 官方 PFT 的保存、全局 ID 映射和恢复流程。

能量统计、DES_USR_VAR(2:31) 和人工 PFR seed 只允许存在于 validation。双精度边界诊断当前作为临时、默认关闭的通用生产钩子保留：仅设置 `MFIX_JKR_DIAGNOSTICS` 时写文件，无固定 case、时间或粒子 ID；完成验收后的验证专用构建迁移列为后续清理项。
