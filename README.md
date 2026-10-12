# 方格棋盘上的 Node Kayles —— 研究工坊

本工坊研究的游戏（正常规则：无棋可走者负）：

> 在 `m × n` 的空棋盘上，双方轮流在空格放一枚棋子；放置时必须保证该格的
> 上、下、左、右四个相邻格中都没有棋子。轮到某人无格可放时，该玩家输。

这个游戏正是**方格图上的 Node Kayles**（也称"独立集游戏"）：棋盘对应图
`P_m □ P_n`（两条路径的笛卡尔积）；双方合力构造方格图的一个极大独立集，
下出最后一手（使独立集极大）的玩家获胜。

**结论总目录：[`RESULTS.md`](RESULTS.md)** —— 工坊里所有的结论（定理、计算值、
新发现、猜想、开放问题），都带有状态标签与证据出处。

## 目录结构

```
RESULTS.md                     结论总目录（从这里开始读）
notes/
  01_game_and_theory.md        游戏规则、Sprague-Grundy 理论、分解原理
  02_symmetry_strategies.md    配对（镜像）引理及其证明与应用
  03_literature.md             文献综述（含出处与核实状态）
  04_computation_log.md        数据如何产生、如何核验
  05_open_problems.md          开放问题与猜想
src/nodekayles/solver.py       参考求解器（纯 Python，易读）
src/cfast/kayles.c             快速求解器（C，位棋盘 + 记忆化）
scripts/run_campaign.sh        重新生成全部数据（幂等、可续跑）
scripts/dawsons_chess.py       1×n 棋盘 = Dawson's chess（octal 0.137）
scripts/analyze.py             序列、周期检测、P-位置、全部校验
scripts/winning_moves.py       小棋盘的首手制胜着（Python）
scripts/pairing_certificates.py 机器可验证的配对策略证书搜索
scripts/play.py                与最优对手交互对弈
tests/test_literature.py       对照已发表数值做验证
tests/cross_check.py           C 求解器 vs Python 求解器，小棋盘全量对照
tests/dumb_solver.py           朴素独立实现（第三方复核）
tests/verify_certificates.py   配对证书的操作性模拟验证
data/                          OEIS 基准数据、原始运行结果、汇总表
```

## 快速上手

```bash
# Python：单块棋盘的 Grundy 值（约 6x8 以内效率尚可）
python3 -c "import sys; sys.path.insert(0,'src'); from nodekayles import grid_value; print(grid_value(5,5))"

# C 求解器：先编译
gcc -O2 -o src/cfast/kayles src/cfast/kayles.c
./src/cfast/kayles board 7 7        # 单块棋盘
./src/cfast/kayles strip 3 30      # 3×n，n=1..30（共享记忆，跨 n 复用）
./src/cfast/kayles wmoves 7 7      # 制胜首着（每轨道取代表格）

# 数据与校验
bash scripts/run_campaign.sh all   # 重算全部数据；已完成的作业自动跳过
python3 scripts/analyze.py         # 序列 + 全部校验（含猜想监测）
python3 scripts/play.py 5 5        # 交互对弈

# 测试
python3 tests/test_literature.py           # 快速档（约 1 分钟）
python3 tests/cross_check.py 6 8           # C vs Python，≤6×8 全量
python3 tests/dumb_solver.py 6 7           # 第三方朴素实现复核
python3 tests/verify_certificates.py       # 配对证书策略模拟
```

## 工坊现状

* 求解器三套：参考（Python）、快速（C）、朴素（独立实现）。三者交叉一致；
  C 版通过"连通分量规范化记忆化 + 对称着法剪枝"比 Python 快约 10³ 倍。
* 数据：`m·n ≤ 约 130` 的全部棋盘，加上 `m ≤ 12` 的各条带家族；各家族的具体
  推进范围见 RESULTS.md（受"状态数按列指数增长"限制）。
* 可机器核验的内容全部核验：OEIS A002187（1001 项）、A316632（17 项）、
  Brown et al. 2020 表 1 的全部可读格，以及三条结构性定理。
* 全部文档为中文；代码注释保留英文。
