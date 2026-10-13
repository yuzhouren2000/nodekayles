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
  06_structural_pairing.md     乘积型配对的完整分类证明与研究路线
  07_fixed_pairing_and_responses.md 任意固定配对反证与 5×8 的应答结构
  08_strip3_reductions.md      三行端口家族、严格剪枝和构造性胜着
  09_strip3_handoff.md         三行下一阶段任务书（可独立接手）
  10_strip3_ports.md           B/E/H 端口族切割/延伸引理、闭合边界、最小反例
  12_two_move_separators.md    两着分割完整分类、小端口零回应与整盘证书
  13_interface_obstruction.md  SG 边界删除表不足以递推的真实可达反例
src/nodekayles/solver.py       参考求解器（纯 Python，易读）
src/cfast/kayles.c             快速求解器（C，位棋盘 + 记忆化）
src/cfast/strip3.c             三行专用求解器（126 位状态 + 严格归约）
scripts/run_campaign.sh        重新生成全部数据（幂等、可续跑）
scripts/dawsons_chess.py       1×n 棋盘 = Dawson's chess（octal 0.137）
scripts/analyze.py             序列、周期检测、P-位置、全部校验
scripts/analyze_strip3.py      3×n 端口族、中行分裂和首着结构检查
scripts/strip3_ports.py        B/E/H 端口族几何、survey/closure/target 调查器
scripts/strip3_two_move.py     用两着分割构造整盘 P 证书、重算辅助端口值
scripts/strip3_interface.py    寻找同边界摘要、接列后胜负不同的局面
scripts/winning_moves.py       小棋盘的首手制胜着（Python）
scripts/pairing_certificates.py 根据结构定理构造配对证书（保留穷举核验）
scripts/fixed_pairing.py        任意固定配对搜索，输出可独立检查的反证
scripts/response_structure.py   首回应调查与固定配对结构提取
scripts/pairing_delay.py        检验两轮应答后能否保证转入固定配对
scripts/play.py                与最优对手交互对弈
tests/test_literature.py       对照已发表数值做验证
tests/cross_check.py           C 求解器 vs Python 求解器，小棋盘全量对照
tests/dumb_solver.py           朴素独立实现（第三方复核）
tests/verify_certificates.py   配对证书的操作性模拟验证
tests/test_strip3_ports.py     B/E/H 切割公式的独立坐标/递推核验
tests/verify_strip3_ports.py   端口族数据文件与最小反例的独立复核
tests/test_strip3_two_move.py  两着分割公式、完整分类与证书负例核验
tests/verify_strip3_two_move.py 独立检查零回应并用原 C 求解器重算 SG 前提
tests/verify_strip3_interface.py 坐标递推与原 C 双重核验接口摘要反例
data/                          OEIS 基准数据、原始运行结果、汇总表
```

## 快速上手

```bash
# Python：单块棋盘的 Grundy 值（约 6x8 以内效率尚可）
python3 -c "import sys; sys.path.insert(0,'src'); from nodekayles import grid_value; print(grid_value(5,5))"

# C 求解器：先编译
gcc -O2 -o src/cfast/kayles src/cfast/kayles.c
./src/cfast/kayles board 7 7        # 单块棋盘
./src/cfast/kayles strip 5 12      # 跨 n 共享记忆
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
python3 tests/test_pairing_classification.py # 分类定理 vs 独立穷举与模拟
python3 tests/test_fixed_pairing.py       # 任意配对搜索 vs 全配对枚举
python3 tests/test_strip3.py              # 3×n 专用求解器的独立核验
python3 tests/test_strip3_ports.py        # B/E/H 切割公式的坐标级核验
python3 tests/verify_strip3_ports.py data/strip3_ports_n12.jsonl   # 端口族数据复核
```

## 当前主线：3×n

三行专用求解器使用 126 位列主序状态：

```bash
make build-strip3
STRIP3_HT_BITS=27 ./src/cfast/strip3 survey 28 > data/strip3_survey.jsonl
make analyze-strip3
make test-strip3
```

调查已精确推进到 **`n=32`（连续）**：`G(3×30)=0`（15 个边行轨道全部非零 +
既有中行排除，见 [notes/10](notes/10_strip3_ports.md) §11），以及本轮的
`f(29)=1、f(31)=2、f(32)=6`（RESULTS.md N20；29/32 另有上边行构造性胜着
`b(29,1)=0`、`b(32,2)=0`）：

```
n :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
G :  2  1  1  0  3  3  2  2  2  3  3  5  2  4  1  3  2  2  1  1  5  0  1  1  3  0  2  5  1  0  2  6
```

`A(n)`（右端中格缺失）与 `D(n)`（两端中格缺失）同时记录，用于验证中行首着
的严格分裂公式；有限前缀的重复只标为观察，不能据此声称终期周期。
完整证明与消融结果见 [notes/08](notes/08_strip3_reductions.md)。已知辅助值直接
给出 `3×32,34,38,40,42,48,52` 的中行制胜证书，也能排除 `3×36` 的全部中行
首着。边行首着之后的端口族 `B,E,H` 已有切割引理与零尾巴延伸引理
（[notes/10](notes/10_strip3_ports.md)），含 `B(10,4)` 的全着法结构、
126 个显式零回应、以及“有限模板词汇表不闭合”的最小反例（T4 爪形）。
新增的[两着分割定理](notes/12_two_move_separators.md)完整分类所有可能的两着
几何分割，并直接给出 `3×4,22,26,30` 的 P 证书。`3×30` 只用宽度不超过 22
的 43 个 SG 前提；原多字 C 求解器约 30 秒即可独立复算其中的 42 个非平凡值，
无需计算任何 B(30,j) 的完整 nimber。运行 `make test-strip3-two-move` 复现。
`3×36` 已有 35/36 轨道回应，边行第 17 列仍未闭合，整盘结果保持未知。

**结构研究的当前边界：** 两着分割提高有限证书的效率，尚未解释空盘零点规律。
[接口摘要审计](notes/13_interface_obstruction.md)进一步发现两个真实可达的
`3×9` 连通 P 局面：末三列及全部 64 项接口删除 SG 表相同，接上一列却分别
得到 SG 0、4。这排除了该摘要作为精确递推状态，未排除更丰富的接口方法。
运行 `make test-strip3-interface` 可独立复算。后续理论目标转向能保留全部
合法后继的接口等价、任意长度替换或空盘的无限胜负族；这些目标尚未完成。

定向求边行后继/整盘值的入口（同一进程内共享记忆，避免逐列重复计算）：

```bash
# 单进程求出某尺寸的全部边行轨道值（轨道数约 n/2）
python3 scripts/strip3_ports.py sweep --n 32 --ht-bits 30 --prealloc \
    --timeout 21600 --output data/strip3_sweep_32.json
# 从 sweep 结果独立重推整盘值 G(3xn)（中行分裂 + 反射 + mex）
python3 tests/verify_strip3_exact.py --n 32 --sweep data/strip3_sweep_32.json
```

状态空间增长律实测（×3.38/2 列）与优化（16 字节打包条目 + 8 线程共享记忆，
约 10× 提速）见 notes/10 §12；union 实测 `U29`=88.9M、`U30`=163.5M、
`U31`=300.7M、`U32`=553.0M，单群 8 线程 3–24 分钟。本机（23GB RAM）表上限
为 2³⁰×16B（7.5×10⁸–9.7×10⁸ 条目），`f(33)` 需要的 `a(32)` union ≈1.4×10⁹
状态，超出本机；n≥34 同样超出（见 RESULTS.md N19/N20）。

原始任务书见 [notes/09](notes/09_strip3_handoff.md)：其阶段 A（L1/L2 证明与
独立测试）、阶段 B（`B(10,4)` 全着法、小族调查、闭合性审计）、阶段 C 的
`3×30` 判定均已完成；`3×36` 是延续目标。

## 工坊现状

* 基础求解器三套：参考（Python）、快速（C）、朴素（独立实现），另有三行专用 C。
  基础三者交叉一致；
  C 版通过"连通分量规范化记忆化 + 对称着法剪枝"比 Python 快约 10³ 倍。
* 数据：`m·n ≤ 约 130` 的全部棋盘，加上 `m ≤ 12` 的各条带家族；各家族的具体
  推进范围见 RESULTS.md（受"状态数按列指数增长"限制）。
* 可机器核验的内容全部核验：OEIS A002187（1001 项）、A316632（17 项）、
  Brown et al. 2020 表 1 的全部可读格，以及三条结构性定理。
* 全部文档为中文；代码注释保留英文。

目前理论优先级见 [notes/13](notes/13_interface_obstruction.md)：检验连接上下文
中的局面等价，防止把有限分割证书当作空盘规律。既有三行分割工具见
[notes/08](notes/08_strip3_reductions.md)。更一般棋盘的配对研究见
[`notes/06_structural_pairing.md`](notes/06_structural_pairing.md)：乘积型固定
配对存在当且仅当“两边均偶数或某边等于 4”，已无须逐尺寸穷举。这不等于
全部 P-位置的分类；`5×8` 等零点仍需研究更一般的配对或动态应答策略。

后续 [notes/07](notes/07_fixed_pairing_and_responses.md) 用完整分支反证排除了
`5×8, 5×12, 7×8` 的任意固定配对，并定位了 `5×8` 的难首着 `(3,2)`（从 1 编号）。
对这一首着，应答两轮后也不能保证转入固定配对；这部分研究尚需带两个中央
接口的动态应答不变量。所有有限负结论附可独立核验的证据文件。
