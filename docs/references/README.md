# JKR validation references

这些 PDF 从只读目录
`/home/yuxing/mfix/P3/JKR/versions/docs/references` 原样复制，用于审查
MFIX 26.1.2 JKR 理论定义和阻尼实现；旧目录未被修改。

| 文献 | 在验证中的用途 | SHA-256 |
|---|---|---|
| Johnson et al. (1971), *Surface energy and the contact of elastic solids* | 经典 JKR 接触半径、表面能和拉脱定义 | `c46dac340ea25c2ca865242c5065bbff782dd820182e20d7d0f8a5db0ae40bba` |
| Marshall (2009), *Discrete-element modeling of particulate aerosol flows* | 恢复系数到法向阻尼参数的拟合来源 | `11fb0e76b22d9eb62ed82effb09565848ba6f7e111391d343c192daffe383719` |
| Chen, Li and Yang (2015), *Sticking/rebound criterion for collisions of small adhesive particles* | 黏附—回弹临界速度与碰撞参数；DOI `10.1016/j.powtec.2015.01.051` | `2c9730e6af0b72084e60fc5a535c4cf3bd60a2def507b8caf48c4b121c93ed7d` |
| Chen et al. (2019), *A fast adhesive discrete element method for random packings of fine particles* | 快速 adhesive DEM 实现背景与交叉核对 | `36bcc9aba8a53b51ca0923d1b6eeb7ac771dffa9c114073ff7c1706025c1e4ab` |

## 使用约束

- 文献公式必须先核对 MFIX26 中 `surface_energy` 与 work of adhesion
  的具体契约，不能把旧 Python 公式直接当真值。
- 当前实现使用 `gamma_eff = sqrt(gamma_i * gamma_j)`，并在经典 JKR
  定义中对应 `W = 2 * gamma_eff`；异质表面能测试会独立验证该组合规则。
- `critical_separation` 是位移控制稳定分支的终点；最大拉力点
  `F = -F_c` 位于另一个接触半径，二者不得混用。
- 未复制 32 MB 的 `MFiX-User_Guide.pdf`，因为当前 conda 环境已提供官方
  MFIX 26.1.2 源码和运行时，且该手册不是 JKR 理论判据的直接来源。
