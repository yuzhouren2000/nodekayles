# 03 —— 文献笔记

本工坊依赖的全部外部结论，附出处与核实程度。
"全文核实" = 下载并阅读了全文（PDF/HTML/OEIS b-file）；
"二手核实" = 只通过引用它的论文确认。

## 游戏本身

术语：本游戏是图 `P_m □ P_n` 上的 **Node Kayles**（也称独立集游戏）。

## 格盘方面的主要论文

**S. Brown, S. Daugherty, E. Fiorini, B. Maldonado, D. Manzano-Ruiz,
S. Rainville, R. Waechter, T. W. H. Wong, "Nimber Sequences of Node-Kayles
Games", Journal of Integer Sequences 23 (2020), Article 20.3.5.**
<https://cs.uwaterloo.ca/journals/JIS/VOL23/Wong/wong24.html>（PDF 全文核实）

* **Theorem 2**：`G(n × 2) = 1`（`n` 奇），`0`（`n` 偶）；用三个辅助"梯子
  变体"做强归纳证明（不是镜像策略）。我们的数据对 `n ≤ 150` 全部验证。
  注意：原文该证明中"（i−1)⊕(n−i) = 1 当且仅当 n 为偶数"一句字面上过宽
  （实际要求 `n ≡ 2 (mod 4)`）；结论本身（奇数 n 时 0 达而 1 不达；偶数 n
  时 0 不达）按原文推导是正确的。
* **Theorem 4**（配对引理）：若图有对合 `φ` 使每个顶点 `v` 的 `φ(v)` 既非
  自身也非其邻居，则 Grundy 值为 0。原文把它同时归功于
  Duchêne–Gravier–Mhalla。
* **Theorem 5**：偶×偶为 0；奇×奇 > 0（占天元 + 180° 旋转）。
* **Theorem 6**：`G(n × 4) = 0` 对一切 `n`。证明用的是"列 j ↔ j+2 (mod 4)"
  的**性质 (P) 不变量**（保持顶点集在这对列置换下封闭的第二手策略）——我们
  核实这**只对恰好 4 列成立**，不能推广到 8 列或 12 列（见 RESULTS.md N8）。
* **表 1**：`1 ≤ n, m ≤ 10` 的 Grundy 值。其中 17 格有底纹（"computed using
  Python"，数值可读），**16 格印刷空白**（作者未计算；代理逐字符核实为白底
  无字）：`(5,8), (5,10), (6,7), (6,9), (7,8), (7,10), (8,9), (9,10)` 及其
  转置。我们已补出其中 10 格（RESULTS.md §2.2），其余 6 格（涉及 9×8、10×7、
  10×9）正在计算。
* **另附**：论文还处理了棱柱图（Corollary 3：`C_n □ P_2` 后手必胜）、链环/
  菱形图等；并明确把 `n × 3` 与一般格图列为开放问题（原文 p.42：
  "such as the n × 3 lattice graph or all lattice graphs in general"）。
  同页还提到 Richard Guy 关于"一切 octal 游戏序列终期周期"的著名猜想。

## 路径情形：Dawson's chess，octal 游戏 0.137

* 在路径上玩 Node Kayles 就是 **Dawson's chess**（T. R. Dawson, 1935），
  octal 游戏 **0.137**：一步可取走 1、2 或 3 枚连续棋子（含所点格及邻格），
  内部着法把整行分裂成两段。识别出处：A. Guignard, E. Sopena,
  "Compound Node–Kayles on paths", Theoret. Comput. Sci. 410 (2009)
  2033–2044（arXiv:0806.3033），全文核实。
* Grundy 序列：**OEIS A002187**（"Sprague-Grundy values for Dawson's Chess
  (octal game .137)"）。OEIS 的表述：*周期 34，例外仅
  `n = 0, 14, 16, 17, 31, 34, 51`*——即 `a(n+34) = a(n)` 对一切 `n ≥ 52` 成立。
  原始文献：R. K. Guy, C. A. B. Smith, "The G-values of various games",
  Proc. Cambridge Philos. Soc. 52 (1956) 514–526；Berlekamp–Conway–Guy,
  *Winning Ways* (1982) pp. 89, 102。
* 本工坊用独立的 O(n²) 区间递推重算：**1001/1001 项与 A002187 完全一致**
  （`scripts/dawsons_chess.py`，基准数据 `data/A002187.txt`）。

## 3 × n 棋盘

* **OEIS A316632**："3 × n 格盘上 Node-Kayles 的 Sprague-Grundy 值"，
  `n = 1..17`：`2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5, 2, 4, 1, 3, 2`；
  唯一零点是 `n = 4`。条目注释记录了 1×n = A002187、2×n = A000035（周期 2）、
  4×n = A000004（恒 0）等识别。未提及任何周期性。
* 我们把它推进到 `n = 27`（还在继续）：新增 P-位置 `n = 22, 26`。
  见 RESULTS.md N1/N7。

## 复杂度

* **T. J. Schaefer, "On the complexity of some two-person perfect-information
  games", J. Computer and System Sciences 16(2) (1978) 185–225**：
  Node Kayles 胜负判定 PSPACE-完全（原文在付费墙后；经
  Demaine–Hearn arXiv:cs/0106019、Grier arXiv:1209.1750、Guignard–Sopena、
  Chandran et al. arXiv:2205.03526 等多方引用交叉确认）。
* **H. L. Bodlaender, D. Kratsch, "Kayles and Nimbers", J. Algorithms 43
  (2002) 106–119**：有界星形数图类上的多项式算法（覆盖路径）。二手核实。
* Y. Kobayashi, "On Structural Parameterizations of Node Kayles",
  arXiv:2003.11775：以点覆盖/模块宽度为参数、邻域多样性核化等 FPT 结果。

## 2026 年的毛毛虫与路径/环幂算法

**M. F. Medeiros-de-Souza, F. E. Hanada, F. Protti, "Impartial Games on
Graphs: Solving Kayles via Dynamic Programming and Periodicity Properties",
Journal of the Brazilian Computer Society 32(1) (2026), 1951–1960.**
[DOI:10.5753/jbcs.2026.7462](https://doi.org/10.5753/jbcs.2026.7462)，
出版日2026-07-31。（已核对出版页面、摘要，未在本轮全文核实PDF。）

* 摘要明确指出：毛毛虫属于区间图，任意毛毛虫上的O(n³)算法已见
  Bodlaender–Kratsch 2002；新工作研究叶子约简和部分毛毛虫族的周期34。
* 另给路径幂、环幂的O(n²)算法。它不是完整三行格盘的求解定理。
* 本工坊独立证明所有三行可达树分量均为毛毛虫，并实现逐手可检查的
  区间递推；2.53秒筛选确认该结构单独没有主线规模收益，见[notes/38](38_caterpillar_components.md)。

## 配对引理的文献形态

* Brown et al. 2020, Theorem 4（原文全文核实）。
* E. Duchêne, S. Gravier, M. Mhalla, "Combinatorial graph games", Ars
  Combin. 90 (2009) 33–44——被引为对合引理的来源（仅经参考文献核实；
  该文本身在网上未能独立找到）。

## 未核实 / 未找到

* 没有找到任何专门研究一般 `m × n` 格盘的论文；该情形在 Brown et al. 2020
  中明确列为开放。
* 我们**没有**找到"后手必胜当且仅当两边都是偶数"这类论断——而且它是**错的**
  （`4 × n` 对奇数 n、`3 × 4`、`6 × 7 = 3` 等都是反例）。请勿引用此类规则。
* Dawson's chess 周期 34 是"谁最先证明"（Guy–Smith 1956 还是后来的计算验证）
  未能从一手资料判定。
* 曾见二次资料称格盘表里 `6 × 7 = 0`：那是 PDF 文本流错位的误读（那个 0 属于
  邻格 `6 × 8`）。我们的三套独立求解器给出 `6 × 7 = 3`，论文原表该格为空白。
  见 RESULTS.md N3。
