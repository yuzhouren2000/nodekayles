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
  14_outer_defect_algebra.md   两行闭合递推、交错尾巴闭式与周期提升
  15_comb_stabilization.md     任意固定核心接交错尾巴的无条件周期定理
  16_clique_column_strips.md   三行特殊列参数递推、完整梯子无界性与周期尾巴
  17_connected_history_language.md 空盘全部连通后继的唯一语言、精确计数与一字节表
  18_boundary_replacement_algebra.md 任意环境的四环—树恒等式与完整端口证书
  19_exact_column_transfer.md 单块生成完整三行图类，接口模 Nim 偏移仍无限
  20_guarded_symmetry.md      固定叶子补偿定理、两手后的密集三行局面 SG 恒为 1
  21_neutral_core_absorption.md 删除中性核心吸收、任意前缀两点尾巴异或恒等式
src/nodekayles/solver.py       参考求解器（纯 Python，易读）
src/cfast/kayles.c             快速求解器（C，位棋盘 + 记忆化）
src/cfast/strip3.c             三行专用求解器（126 位状态 + 严格归约）
src/cfast/strip3_language.c    全可达连通家族的完美编号求解器（每形状一字节）
scripts/run_campaign.sh        重新生成全部数据（幂等、可续跑）
scripts/dawsons_chess.py       1×n 棋盘 = Dawson's chess（octal 0.137）
scripts/analyze.py             序列、周期检测、P-位置、全部校验
scripts/analyze_strip3.py      3×n 端口族、中行分裂和首着结构检查
scripts/strip3_ports.py        B/E/H 端口族几何、survey/closure/target 调查器
scripts/strip3_two_move.py     用两着分割构造整盘 P 证书、重算辅助端口值
scripts/strip3_interface.py    寻找同边界摘要、接列后胜负不同的局面
scripts/strip3_defect_algebra.py 两行区间代数、任意核心交错尾巴的周期求值
scripts/periodic_clique_strip.py 任意三行有限缺陷递推及混合周期尾巴求值
scripts/strip3_history_language.py 全可达连通局面的规范编码、计数、编号与 SG
scripts/guarded_symmetry.py    受保护固定集的精确 SG 证书与密集三行后继族
scripts/neutral_core.py        核心删除中性与配对外围证书、局部尾巴直接消去
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
tests/verify_comb_stabilization.py 穷尽与核心大小无关的四值块稳定性引理
tests/test_clique_column_strips.py 任意三行缺陷代数、端点表与真实代表着法核验
tests/verify_clique_strip_period.py 三行混合尾巴周期窗口的原图独立复算
tests/test_strip3_history_language.py 全部历史/组件、双向编号及字节表独立核验
tests/test_guarded_symmetry.py 固定集奇偶公式、完整着法闭合与真实历史独立核验
tests/verify_neutral_core.py   独立检查核心全部 mex/删除状态及外围配对条件
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
运行 `make test-strip3-interface` 可独立复算。

新增[外行缺口代数](notes/14_outer_defect_algebra.md)给出正面结果：第三行消失后，
任意两行诱导局面的全部着法落入 `O(n³)` 时间的精确区间递推。交替缺口尾巴
有闭式周期 4；以单条边连接任意固定头部时，对全部保根后继核验连续六处
周期比较，便可严格归纳到无限长度。6 个可达三行实例的全部 50,407 项不同
SG 前提已由原 C 求解器重算。特别地，在空 `3×(4k+6)` 上下出上行
`(1,5),(1,7),…,(1,4k+5)`，所得连通局面为 P，对一切 `k≥0` 成立。
周期证书还将任意长尾巴的制胜着搜索缩到有限代表；运行
`python3 scripts/strip3_defect_algebra.py --head 4 --tail 1000000` 得到 SG 2
与制胜着 `(3,1000004)`，只需 26 个候选。`make test-strip3-defect` 独立复现。
这些是缺口局面的无限规律，空三行棋盘的完整分类仍未得到。

进一步的[无条件稳定定理](notes/15_comb_stabilization.md)证明：任意 `k` 顶点
核心以单边接交错尾巴，全部 SG 都小于 `4(k+2)`，并且最迟从
`192k²+984k+3` 以内的某个起点进入周期 4。四值块有限引理加上历史集合
饱和、二值 mex 映射给出证明，适用于所有有限核心，核验窗口必然能找到。
因此 `RootedComb.stabilize()` 无须猜测尾巴截断长度；
`make test-comb-stabilization` 核验通用引理及原有完整选项谱。

进一步的[三行特殊列代数](notes/16_clique_column_strips.md)把任意三行诱导图
压成区间、端点删除和 h 个特殊列的三态选择，得到 `O(3^h n³)` 精确递推；
h 是原来同时含上下外格的列数。纯外行历史只要覆盖每列及邻列，就有 h=0，
第三行可反复保留，SG 仍能多项式求出。混合周期历史“上、上、下、下”，列距
`2,3,2,3`，给出第三行仍存的周期10尾巴及连接三行头部的无限 P 族。
另一方面，完整两行梯子删一个角格的 SG **恰为 n**，证明周期几何和固定核心
不保证 SG 有界或最终周期。14,989 个周期前提由原 C 独立复算；任意 `3×5`
诱导图及接口状态由朴素原图递推核验。运行 `make test-clique-strips`。
空盘 h=n，仍未得到完整空三行棋盘的分类或周期。

[连通历史语言定理](notes/17_connected_history_language.md)进一步去掉形状
限制：空三行盘的**每个可达连通分量**都有唯一“正则外行词 + 两端有限删除”编码。
同行落子间距至少 2，异行至少 3；全跨连通后继数严格满足三阶递推，全部连通
掩码数也从第 7 项起满足同一递推，指数底数 `ρ=1.839286755…`。它解释了旧
求解器每两列约 3.38 倍的增长，也把 SG 表改成无哈希、每几何一字节：完整跨度
至 33 的表容量为 9.21 GB，已落入本机内存范围，**但 f(33) 仍未求得**。
新 C 版单线程查询较慢，测试已复现到 22；运行 `make test-strip3-language`。
内部落子列都是割点，所有环局限在未落子段；厚块之间共享割点的游戏接口仍待解决。

最新的[连接代数](notes/18_boundary_replacement_algebra.md)保留全部内部后继，
证明真实三行四环图可在任意指定双端口环境下换成树。其
[逐列运算](notes/19_exact_column_transfer.md)进一步用一个六点 `3×2` 生成元
表示任意三行诱导图，完整空盘成为同一列算子的轨道。同时证明全部诱导图的
连接接口即使模去任意独立 Nim 堆也仍无限：缺角梯子连接延长实行整数加法，
不能由固定异或偏移补偿。实现仍指数；这一障碍未排除空盘实际可达家族的
特殊结构。运行 `make test-kayles-transfer`；没有新空盘数值或周期结论。

新的[受保护固定集定理](notes/20_guarded_symmetry.md)给出精确消去规律：对合
固定集中的每个中心若带奇数个全图固定叶子，其余点成非相邻镜像对，则
SG 等于中心数与固定孤立点数之和的奇偶。镜像回应删除中心时，幸存叶子的
奇偶恰好补偿；外围的三行厚块、环和跨侧边均不需要求值。任意奇数 `n≥3`
的空三行盘，在上行中央左右各一格落子后，所得连通局面 SG 恒为 1，中央
中、下两格都是明确制胜着；`n≥5` 时仍保留 `3n−7` 格。它还给出所有
`n≥8` 的四手 P 历史，但不声称能强迫对手采用这些历史，也未判定空盘。
运行 `make test-guarded-symmetry`；证明及独立核验见 notes/20。

**删除中性核心与局部改写：** [notes/21](notes/21_neutral_core_absorption.md)
将固定叶子补偿推广成遗传的端口删除等值条件。满足条件的核心可以吸收任意
兼容配对外围，精确保持核心 SG；配对无需图自同构，核心值可取任意 nimber
（一般图构造，未声称全部扩图三行可达）。另一条纯结构规则证明：已有叶子的
顶点再接两点尾巴，SG 恰好异或 1。因此任意三行前缀 U 都有
`g(U,4,6,4,4)=g(U,4,6) xor1`；真实三手后继只删除 12 格也能应用，但剩余
厚核心尚未解出。运行 `make test-neutral-core`；核心和外围条件均可独立验算。

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

目前完整可达几何的结构与下一道边界见
[notes/17](notes/17_connected_history_language.md)：连通历史已有唯一语言、
精确计数和完美编号；SG 的无限规律、厚块的共享割点接口仍待研究。
任意上下文替换及完整列运算见 [notes/18](notes/18_boundary_replacement_algebra.md)、
[notes/19](notes/19_exact_column_transfer.md)；后者也划出通用有限接口压缩的边界。
任意对称外围的固定叶子补偿及密集三行精确值见
[notes/20](notes/20_guarded_symmetry.md)；空盘如何进入这些已解结构仍待研究。
更一般的核心吸收与任意前缀尾巴恒等式见
[notes/21](notes/21_neutral_core_absorption.md)。
任意三行有限缺陷与混合周期尾巴见 [notes/16](notes/16_clique_column_strips.md)。
任意固定核心接交错尾巴的无条件稳定性
见 [notes/15](notes/15_comb_stabilization.md)。
接口设计的反例基准见 [notes/13](notes/13_interface_obstruction.md)。既有三行分割工具见
[notes/08](notes/08_strip3_reductions.md)。更一般棋盘的配对研究见
[`notes/06_structural_pairing.md`](notes/06_structural_pairing.md)：乘积型固定
配对存在当且仅当“两边均偶数或某边等于 4”，已无须逐尺寸穷举。这不等于
全部 P-位置的分类；`5×8` 等零点仍需研究更一般的配对或动态应答策略。

后续 [notes/07](notes/07_fixed_pairing_and_responses.md) 用完整分支反证排除了
`5×8, 5×12, 7×8` 的任意固定配对，并定位了 `5×8` 的难首着 `(3,2)`（从 1 编号）。
对这一首着，应答两轮后也不能保证转入固定配对；这部分研究尚需带两个中央
接口的动态应答不变量。所有有限负结论附可独立核验的证据文件。
