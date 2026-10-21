# MFIX 26.1.2 球形颗粒 JKR

本项目以 MFIX 26.1.2 官方源码为基底，实现球形颗粒的 JKR 黏附接触、滑移、
滚动和扭转，并维护邻居重建、并行迁移及续算中的接触历史。
当前本地源码基线为 `9b5d639`，对应 GitHub 源码发布
[`e72b226`](https://github.com/yuxing081/MfixVer26/commit/e72b22660868c2c4ac9d9ed16013f87659fadf35)。
`code/` 保存相对官方版本的同名覆盖文件，不是完整的 MFIX 安装包。

## 相对官方 MFIX 的主要修改

| 修改 | 当前实现 | 主要代码入口 |
|---|---|---|
| 球形 JKR 接触与阻力 | 颗粒对及颗粒–墙面的法向黏附、切向滑移、弹性滚动、滚动阻尼和扭转；支持相、墙面和单颗粒表面能，以及实际/模拟杨氏模量对应的降刚度参数映射 | [颗粒对](code/calc_force_dem.f)、[墙面](code/calc_collision_wall_mod.f)、[共享内核](code/usr_mod.f)、[参数检查](code/check_solids_dem.f) |
| 接触历史方向 | 在法向转动和共同自转时运输切向历史 `PFT` 与滚动历史 `PFR`，再更新增量及限幅；`PFTW` 保存相对扭转角；移除近零历史因舍入残量触发的误停机判据 | [历史运输](code/usr_mod.f)、[历史变量](code/discretelement_mod.f) |
| MPI 迁移与邻居重建 | 沿用 MFIX 的全局颗粒编号和接触槽，补齐任一端点迁移、角点转发及历史恢复；发生所有权迁移时触发邻居重建，使历史在下一次力计算前可用 | [迁移通信](code/mpi_funs_des_mod.f)、[打包](code/mpi_pack_des_mod.f)、[解包](code/mpi_unpack_des_mod.f)、[邻居管理](code/neighbour.f) |
| 壁面与邻居几何 | 在 JKR 路径处理共面 STL 接缝的重复接触和安全历史转移，并修正半网格邻居剪枝阈值 | [墙面接触](code/calc_collision_wall_mod.f)、[网格搜索](code/desgrid_mod.f) |
| 重启与停止保存 | 保存 JKR 接触状态、滑移/滚动/扭转历史、旧法向和历史时间；当前 JKR 重启格式为 2.2；纯 DEM 收到停止请求时写出成对重启文件 | [重启读取](code/read_res0_des.f)、[重启写出](code/write_res0_des.f)、[停止处理](code/mfix.f) |
| 力链输出与 GUI 输入 | 按输出区域选择接触状态、接触半径、有效表面能、切向位移、滚动位移和弹性扭转角六个字段；归一化使用副本，避免多区域输出覆盖原始力幅值；注册项目关键字 | [VTP 输出](code/vtp_mod.f)、[输出参数](code/output_control.inc)、[GUI 关键字](code/usr_init_namelist.f) |

JKR 由 `USE_JKR=.TRUE.` 启用，使用普通球形 Hertz 接触路径
（`DES_COLL_MODEL='HERTZIAN'`）。不支持非球形 SuperDEM/GSP、LSD，
也不能与原生 `USE_COHESION` 混用。JKR 专用物理和历史扩展受开关控制；
多区域输出与纯 DEM 停止保存修复也适用于关闭 JKR 的情况。

已有验证记录覆盖接触响应、历史方向、壁面几何、特定 MPI 分解与迁移工况、
同版本续算和短时千颗粒剪切。上述结果不代表任意进程数、改变分区后的续算，
或长时间团聚体破碎与宏观流变已得到完整验证；本次发布没有重新运行求解器。

## 编译与运行

在已安装 MFIX 26.1.2 的 WSL/Linux 环境中，从本仓库的 `code/` 目录构建：

```sh
conda activate mfix-26.1.2
cd /path/to/MfixVer26/code
build_mfixsolver -j        # 串行
build_mfixsolver --dmp -j  # MPI/DMP
```

在算例目录准备 `mfix.dat` 和所需颗粒输入，调用本次构建的入口：

```sh
cd /path/to/case
/path/to/MfixVer26/code/mfixsolver
# DMP 示例：输入中的分区乘积应与进程数一致
mpiexec --bind-to none -np 8 /path/to/MfixVer26/code/mfixsolver_dmp
```

构建通过同名文件覆盖机制使用 `code/` 的源码，保持该目录的扁平布局。
不修改官方安装源码，不使用旧构建目录中的残留求解器。需要保存接触历史的
续算使用 `RESTART_1`，并同时保留同次运行的 `RUN_NAME.RES` 和
`RUN_NAME_DES.RES`；当前 JKR 2.2 读取路径不兼容旧 2.0/2.1 历史文件。

## 文档入口

以下是本地完整工作区的文档导航。本次 GitHub 更新仅发布源码和本 README；
其他文档保留远程已有历史版本，尚未同步的本地文档以路径列出。
`validation/` 是独立仓库，直接克隆本仓库不会包含其验证记录或算例。

| 需要了解的内容 | 文档 |
|---|---|
| 当前能力、已验证范围、已知限制 | `docs/STATUS.md` |
| 公式、参数、降刚度、滚动、扭转和历史方向 | `docs/CONTACT_MODEL.md` |
| 接触身份、MPI 迁移、墙面历史和重启格式 | `docs/CONTACT_HISTORY.md` |
| 编译、时间步、停止续算、力链字段和诊断 | `docs/RUNNING_AND_OUTPUT.md` |
| 原生缺陷的原始证据与当时处理 | [MFIX 原生实现不足](docs/MFIX_UPSTREAM_LIMITATIONS.md) |
| 文献及外部实现参考（本地索引较远程新） | [参考资料](docs/references/README.md) |
| 迁移过程、旧版本结论、原文恢复方法 | `docs/archive/CHANGE_RECORD.md` |

原生实现不足文档保留历史记录，其中的旧状态以本页及本地 `docs/STATUS.md` 为补充。
详细验证方法与容差见本地 `validation/docs/VALIDATION.md`，阶段进度见
`validation/docs/TEST_PROGRESS.md`。

## 目录与职责

下表描述本地完整工作区；除 `code/` 和远程已有文档外，所列内容不一定随本仓库发布。

| 目录 | 用途 |
|---|---|
| `code/` | 相对官方版本修改的生产 Fortran 文件，按原生同名覆盖机制保持扁平布局 |
| `model/` | 官方只读源码快照，不纳入生产仓库版本控制 |
| `docs/` | 当前技术说明、资料索引和历史记录 |
| `tests/force_chain/` | 生产输出接口回归 |
| `validation/` | 独立 Git 仓库，维护算例、参考公式、验证入口及结果 |
| `UsrVer26/` | 独立 Git 仓库，维护各用户定制版本及同步记录 |
| `FlowVer26/` | 流动算例和后处理工作区 |

只在对应目录维护其职责内容。生产代码不放算例专用历史种子或能量统计；
现有双精度诊断由 `MFIX_JKR_DIAGNOSTICS` 显式启用，默认关闭。

## 基准与运行约定

- 官方环境：`/home/yuxing/miniforge3/envs/mfix-26.1.2`，在 WSL Ubuntu-22.04 中运行。
- 官方源码：`/home/yuxing/miniforge3/envs/mfix-26.1.2/share/mfix/src`，不在该目录构建或修改。
- 旧 JKR 来源：MFIX 23.4.1，提交 `ea55027498a4b939a2982cbb77d137201b034ffc`；只迁移经核对的模型与历史语义。
- 当前 JKR 重启格式为 2.2；不能用改写文件头的方式把旧历史变成兼容状态。

完整编译和验证说明维护在本地 `docs/RUNNING_AND_OUTPUT.md`。开始新工作前
分别检查生产、验证和用户版本仓库的 Git 状态，保留已有改动。
