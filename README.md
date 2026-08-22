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
- validation/cases/：按 particle_wall、particle_pair、wall_geometry、dmp 和 restart 分类；输入由集中配置在运行时生成。
- validation/jkr_tools/：配置、case 生成、独立理论参考、I/O、比较和绘图工具。
- validation/tests/：不依赖 solver 的单元测试。
- validation/results/：每次运行的独立时间戳目录；大型输出默认不进入 Git。
- docs/VALIDATION.md：正式验证设计、执行结果和限制文档。
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

    cd /home/yuxing/mfix/MfixVer26
    ./validation/run_validation.sh unit
    ./validation/run_validation.sh restart
    ./validation/run_validation.sh all

当前已接入 `unit`、官方 stock 对照的串行 `hertz` 与串行 `restart`；其余 JKR 物理、速度、几何和 DMP 套件按 `docs/VALIDATION.md` 的矩阵逐阶段接入，未运行项目不会标记为通过。

restart 套件对两颗粒和斜碰墙面分别运行：

1. continuous；
2. split stage + RESTART_1；
3. split stage + RESTART_2 敏感性负对照。

RESTART_1 必须与 continuous 的双精度 `particle_output.csv` 逐字节一致。RESTART_2 是否不同单独记录，不作为正向通过条件。结果写入 `validation/results/runs/<UTC时间戳>-<pid>/`，永不复用固定目录。
    ./validation/run_validation.sh hertz

完整测试顺序保持为原生 Hertz 保护、串行墙面/颗粒对、速度扫描、壁面几何、2-rank DMP、RESTART_1 和小颗粒群集成测试。正式结果写入 `docs/VALIDATION.md`，活动进度写入 `docs/TEST_PROGRESS.md`。

## 生产代码与验证代码边界

生产 code/ 保留 JKR 力学、PFT/PFR、JKR contact state、逐颗粒表面能、邻居/DMP/restart 生命周期所必需的逻辑。新增历史仿照 MFIX 官方 PFT 的保存、全局 ID 映射和恢复流程。

能量统计、DES_USR_VAR(2:31)、人工 PFR seed、ghost/slot/force 导出和双精度边界快照只允许存在于 validation 的具体 case、jkr_tools 或显式验证 patch 中，不得无条件进入生产力计算。
