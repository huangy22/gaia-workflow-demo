超导比热计算
基于 17 篇文献的工作流分析
May 11, 2026

1

引言

超导体（superconductor）的比热（specific heat）是刻画其低能准粒子（quasiparticle）谱
与配对对称性（pairing symmetry）的关键热力学观测量之一。与电输运（transport）或
磁响应（magnetic response）相比，比热直接关联体系的态密度（density of states, DOS）
与熵（entropy）
，因此对能隙结构（energy gap structure）
、节点（node）/全能隙（full
gap）特征、以及外场（magnetic field）与杂质（impurity）引入的低能态尤为敏感。在
非常规超导（unconventional superconductivity）研究中，低温比热的幂律（power law）
或指数（exponential）行为常被用于区分结点型与全能隙型配对，并为多带（multiband）
或各向异性（anisotropic）超导提供约束。
另一方面，随材料体系从传统 BCS（Bardeen–Cooper–Schrieffer）到强耦合（strong
coupling）、多轨道（multi-orbital）与拓扑超导（topological superconductivity）拓展，
比热理论计算的输入与实现路径也日益多样：既可能从微观哈密顿量（microscopic
Hamiltonian）出发，经由能谱（spectrum）求得自由能（free energy），也可能直接
以实验拟合或数值模拟给出的 DOS、能隙函数（gap function）为输入构建热力学
势（thermodynamic potential）。此外，在有限温度（finite temperature）、有限外磁
场或涡旋态（vortex state）背景下，如何在保证自洽（self-consistency）与可控近似
（controlled approximation）的前提下稳定获得比热，涉及到热力学求导（thermodynamic
differentiation）的数值敏感性、低温极限（low-temperature limit）的处理以及特征量
（characteristic quantities）提取等一系列技术选择。围绕上述问题，相关文献已形成多
条常用路线，本文汇总并梳理了共 17 篇代表性工作，旨在形成可复用、可比对的“超
导比热计算”工作流视角。
本文采用工作流方式对超导比热计算进行综述，核心流程可概括为四个阶段：（i）
谱/DOS/自由能输入：根据研究目标与可得信息，选择以准粒子能谱（quasiparticle
spectrum）、态密度（DOS）或自由能（free energy）作为起点输入；（ii）构建热力学
势/熵：由配分函数（partition function）或有效作用量（effective action）得到热力学
势（thermodynamic potential，例如自由能）并进一步得到熵（entropy）；（iii）热力学
求导得比热：在保持适当自变量（如体积、粒子数或化学势（chemical potential））不
变的条件下进行求导，获得定容比热（heat capacity at constant volume, CV ）或定压比
热（heat capacity at constant pressure, CP ）
，并区分电子贡献（electronic contribution）
与声子贡献（phonon contribution）
；
（iv）极限/数值评估与特征量：在低温、临界温度
（critical temperature, Tc ）附近或特定外场极限下评估结果，提取比热跃变（specific-heat
jump, ∆C）
、Sommerfeld 系数（Sommerfeld coeﬀicient, γ）
、残余线性项（residual linear
term, γ0 ）等特征量，并讨论数值稳定性与误差来源。为便于读者在不同输入条件与目
标输出之间快速定位合适路线，本文给出对应的工作流决策树（decision tree）示意，如
图 1 所示。

1

Figure 1: 工作流决策树
在此框架下，本文将不同方法学分支纳入统一叙述：从 BCS 平均场（mean-field）
与 Bogoliubov–de Gennes（BdG）谱出发的直接求和/积分路线；以格林函数（Green’s
function）或准经典理论（quasiclassical theory，例如 Eilenberger/Larkin–Ovchinnikov
方程）获取 DOS 并重建熵与比热的路线；以及以自由能泛函（free-energy functional）
2

或有效场论（effective field theory）为中心、通过温度导数得到热容的路线。通过对这
些路线在输入需求、近似层级、数值代价与可解释性方面的对比，本文旨在为不同材料
体系与数据条件下的比热计算提供可操作的选择依据与交叉验证思路。

谱/DOS/自由能输入

2

本阶段目标是为后续热力学量（熵、比热、自由能差等）的计算提供自洽的谱学/热力
学输入：准粒子能谱 Ek 、态密度（density of states, DOS）N (E)，或热力学势/自由能
F (Ω)。在弱耦合 BCS 类框架中，常以已知或自洽求得的 N (E) 作为起点，结合费米–狄
拉克分布构造超导态熵与比热的可计算表达式；当考虑各向异性（含点/线节点）、多
带、杂质与磁场时，差异主要体现在 Ek 与 N (E) 的低能行为及其角向平均方式上。

2.1

由 DOS 构造热力学输入：熵与比热的统一起点

参考文献提取（[6]）给出了以 N (E) 为核心输入的热力学计算流程。设
• E ≥ 0：准粒子能量变量；
• N (E)：单位体积、单位能量的准粒子态密度（超导态）；
• fFD (E, T )：费米–狄拉克分布；
• γST ：正常态线性比热系数标度（normal-state Sommerfeld coeﬀicient，用于归一
化）；
• Tc ：超导临界温度。
费米–狄拉克分布为
fFD (E, T ) =

1
eE/(kB T ) + 1

,

(1)

其中 kB 为玻尔兹曼常数。
（1）以 N (E) 计算超导态熵 Ss (T ) 对具有自旋简并等因素的准粒子系统，其超导态
熵可写成（[6]）
Z ∞
h

i
Ss (T ) = −4
dE N (E) fFD ln fFD + 1 − fFD ln 1 − fFD ,
(2)
0

其中系数 4 吸收了自旋与粒子–空穴等计数约定；在具体实现中应与所采用的 N (E) 规
范保持一致（例如 N (E) 是否已包含自旋简并）
。
（2）由熵得到超导态比热 Cs (T ) 并归一化 比热通过


∂Ss
Cs (T ) = T
∂T

(3)

获得，并常用 γST 进行归一化以比较不同材料或不同相的热力学特征：
Cs (T )
γST T

或

3

Cs (T )
.
γST

(4)

在 T = Tc 处可读出比热跃迁幅度（示意写法）
∆C
Cs (Tc− ) − Cn (Tc+ )
=
,
C Tc
Cn (Tc+ )

(5)

其中 Cn 为正常态比热（常近似 Cn = γST T ，但应与实验/模型一致）。

2.2

各向异性能隙与低能 DOS：低温渐近作为一致性检验

对于各向异性超导（含点节点/线节点）
，能隙角向结构会决定低能 N (E) 的幂律行为，从
而决定低温比热的幂律。文献提取指出：点节点各向异性能隙在低能产生相应的 DOS
行为，导致低温比热满足
Cs (T ) ∝ T 2 ,
(T ≪ Tc ),
(6)
该关系常被用作检验所输入/自洽得到的 N (E) 是否与设定的节点结构相一致。

2.3

自洽能隙求解对 DOS 的支撑（弱耦合、相依赖角因子）

在许多实现中，N (E) 并非完全经验输入，而由自洽能隙 ∆(T )（及其角向形式）确定。
文献提取（[6]）提到针对两种相（A/B 相）输入角向形式因子 f ，并在自洽方程中施加
能量/频率截断 E0 。一种常见的抽象写法是
Z E0
∆(T ) = λ

D
E
,
dξ K ξ, ∆(T )f (k̂), T
k̂

0

(7)

其中
• λ：有效耦合常数（弱耦合参数）；
• ξ：相对费米能的能量变量；
• E0 ：能量（或 Matsubara 频率）截断；
• f (k̂)：能隙角向形式因子（由相 A/B 给定）；
• ⟨· · · ⟩k̂ ：对费米面方向的平均；
• K：由所选理论（BCS/Eilenberger 等）决定的核函数。
在得到 ∆(T ) 与 f (k̂) 后，Ek 与 N (E) 可按所选模型计算并回代到式(18)–(24)中形成闭
环。

4

推荐工具
常用工具：
• 数值积分器（自编或通用库）
：用于计算式(18)的能量积分与角向平均；建议
支持自适应积分与高精度浮点（用途：抑制低温 fFD 的数值病态）。
• 非线性方程求解器（牛顿法/割线法/Brent 法）
：用于求解自洽能隙方程（用
途：得到 ∆(T ) 并生成一致的 N (E)）。
• 数据处理与绘图软件（未在文献中报告具体软件）：用于生成 Ss (T )、Cs (T )
曲线与在 Tc 处读取 ∆C/C（用途：参数扫描与可视化对比）。
文献记录的工具信息（[6]）：
• 数值积分与数值求解：用于 ∆(T )、Ss (T ) 与 Cs (T ) 曲线计算；具体软件未
报告。
关键参数
关键参数与调优建议：
• 能量积分上限：理论上为 [0, ∞)；数值实现建议取 Emax 使得 fFD (Emax , T ) ≪
1（例如对最低温 Tmin 令 Emax 20 kB Tmin ），以避免截断误差进入 Ss 。
• 自洽方程截断 E0 ：需与弱耦合模型一致并做收敛性检查；建议对 E0 做倍增
扫描（例如 E0 → 2E0 ）验证 ∆(T ) 与 Cs (T ) 变化可忽略。
• 能量网格：低温下 N (E) 与 fFD 在 E ∼ kB T 附近最敏感；建议在 E10 kB T
区间加密网格（对数网格或自适应网格），并检查 Cs (T ) 的网格收敛。
• 角向平均（各向异性/节点）：对 f (k̂) 的节点附近需加密采样；建议使用分
层采样或重要性采样以稳定得到低能 N (E) 幂律，从而验证式(6)。
• 归一化参数 γST ：用于与实验或不同相对比；建议明确其来源（实验正常态
拟合或模型计算），并在报告中给出单位与采用的规范（是否含自旋简并）。

5

常见陷阱
常见错误与解决方案：
• 忽略 N (E) 的数值噪声/非物理解：低温 Ss 对 N (E ≈ 0) 极敏感，微小噪声
会放大到 Cs 的振荡。解决：对 N (E) 做物理约束（非负、平滑）与网格收
敛；必要时用保形插值/轻度正则化并报告影响。
• 能量积分截断不足：Emax 太小会低估熵并导致 Cs 系统性偏差。解决：按
Emax 20 kB Tmin 设定并做 Emax 收敛测试。
• 自洽截断 E0 与热力学积分尺度不一致：用很小的 E0 得到的 ∆(T ) 可能
导致与实验不符的 Tc 与 ∆C。解决：对 E0 做扫描并同时检查 Tc 、∆(0)、
∆C/C 的稳定性。
• 简并因子（式(18)前因子）与 N (E) 定义不匹配：会造成整体幅度错误。解
决：在实现中固定约定：N (E) 是否为单自旋 DOS、是否包含粒子–空穴计
数，并据此调整前因子。

构建热力学势/熵

3

本阶段目标是从统计物理出发，将微观输入（如准粒子态密度 N (E)、能谱 Ek 、配分函
数 Z/Ω 或费米分布）转换为可与实验直接对比的热力学量，核心产物为超导态熵 S(T )
（以及由此导出的比热 C(T )、内部能 U (T )、自由能 F (T ) 等）。在提取文献中，方法路
径以 [6] 的 “由 N (E) 与费米分布构造 Ss (T )，再求 Cs (T ) 并归一化到 γST ” 为代表。

3.1

由准粒子 DOS 构建熵与比热（以 [6] 为代表）

在弱耦合准粒子图像下，给定超导态准粒子态密度（density of states, DOS）N (E)，超
导态熵可由费米-狄拉克分布的二元熵积分得到。令
fF D (E, T ) =

1
eE/(kB T ) + 1

,

(8)

其中 kB 为玻尔兹曼常数，T 为温度。文献 [6] 采用的超导态熵（每单位体积或每摩尔
的归一化由 N (E) 的定义决定）为
Z ∞
h

i
Ss (T ) = −4
dE N (E) fF D ln fF D + 1 − fF D ln 1 − fF D .
(9)
0

式中：
• Ss (T )：超导态熵；
• N (E)：超导态准粒子 DOS（能量 E 相对费米能的激发能，取 E ≥ 0）；
• 前因子 4：通常来自自旋简并（×2）以及正负能量对称性下从 (−∞, ∞) 折半到
(0, ∞)（×2）；具体约定需与 N (E) 的定义一致。

6

得到 Ss (T ) 后，比热由热力学关系给出：
Cs (T ) = T

∂Ss (T )
.
∂T

(10)

在与实验对比时，[6] 进一步使用正常态的线性比热系数（Sommerfeld 系数）γST 进行
标度或归一化，典型输出包括 Cs (T )/γST 以及在临界温度 Tc 附近的比热跃迁
∆C
C Tc

或

∆C
,
γST Tc

(11)

其中 ∆C = Cs (Tc− ) − Cn (Tc+ )，Cn (T ) ≈ γST T 。

3.2

低温渐近与节点结构的约束

若能隙具有点节点并导致低能 DOS 呈幂律行为，[6] 指出可推得低温比热的幂律渐近。
例如当低能态密度满足
N (E) ∝ E (E → 0),
(12)
则由式 (9)–(10) 可得到
Cs (T ) ∝ T 2

(T ≪ Tc ),

(13)

从而用 Cs (T ) 的幂律来反推能隙的节点各向异性特征。

3.3

与能隙自洽求解的耦合（支撑 N (E) 的生成）

上述热力学计算依赖输入 N (E)。在 [6] 中，N (E) 由弱耦合自洽能隙方程支撑，并考虑
两种相（A/B）对应的角向形式因子（记为 f (k̂)）。概念上可写作
∆(k̂, T ) = ∆(T ) f (k̂),

(14)

并在能量（或 Matsubara 频率）求和/积分中施加截断
E0 以表征有效相互作用能窗。
q
该自洽部分决定了准粒子谱 Ek =
得到 Ss (T ) 与 Cs (T ) 曲线。

ξk2 + |∆(k̂, T )|2 ，进而决定 N (E)，最终通过式 (9)

推荐工具
常用工具：
• 数值积分工具：用于计算式 (9) 的能量积分（E ∈ [0, ∞)）以及温度扫描；需
要支持自适应积分或高密度网格积分以处理 fF D 在低温下的陡峭变化（用
途：得到 Ss (T )、Cs (T )）。
• 数值求解器（非线性方程/自洽迭代）
：用于弱耦合自洽求解 ∆(T ) 并生成与
相（A/B）角向形式因子 f 一致的 N (E)（用途：建立 N (E) → Ss → Cs 的
闭环计算）。
•（文献 [6] 未报告具体软件）
建议软件：
Python+NumPy/SciPy
（quad/quadpack、
root）、MATLAB（integral、fsolve）、或 Julia（QuadGK、NLsolve）（用
途：快速复现实用的积分与自洽迭代管线）。

7

关键参数
关键参数推荐值与调优建议：
• 能量上限 Emax （用于以有限区间近似 [0, ∞)）：建议取 Emax 20 kB Tmax ，以
保证费米尾部对式 (9) 的贡献可忽略；若 N (E) 在高能仍缓慢变化，可适当
增大至 50 kB Tmax 并检查收敛。
• 能量网格与自适应策略：低温时 fF D (E, T ) 在 E ∼ kB T 附近变化最剧烈，
建议对 E ∈ [0, 10kB T ] 加密采样或使用自适应积分；以 Ss (T ) 在网格加密后
相对变化 < 10−4 作为收敛判据。
• 温度微分求 Cs (T )：由式 (10) 数值求导时，建议使用中心差分
∂S
S(T + δT ) − S(T − δT )
≈
,
∂T
2δT
并取 δT /T ∼ 10−3 –10−2 ；同时对不同 δT 做稳定性检查，避免导数噪声导
致的虚假振荡。
• 自洽截断能 E0 ：[6] 提到在自洽方程中施加能量/频率截断 E0 ；建议将 E0
作为模型超参数，与 Tc 、∆(0) 的拟合联合标定，并检查 E0 改变时低温幂
律（如 T 2 ）是否稳健。
• 归一化系数 γST ：用于输出 Cs /γST 并读取 Tc 处跃迁；建议优先采用同一样
品、同一拟合窗口得到的正常态 γ，避免因背景扣除不一致造成 ∆C 偏差。
常见陷阱
常见错误与解决方案：
• 前因子与 N (E) 约定不一致：式 (9) 的 “4” 依赖于自旋简并与能量对称性的
计数。解决：明确 N (E) 是否已包含自旋因子、是否定义在全能轴；用正常
态极限（高温或 ∆ → 0）下 Cn = γT 的一致性来校验整体因子。
• 将 [0, ∞) 积分粗暴截断导致 S(T ) 偏小：若 Emax 过低，会漏掉中高能尾部
贡献。解决：用 Emax 扫描做收敛测试，直至 Ss (T ) 与 Cs (T ) 变化可忽略。
• 数值求导放大噪声，导致 Cs (T ) 出现非物理振荡：解决：对 Ss (T ) 先做平
滑/样条拟合再求导，或直接对积分表达式做解析求导（将 ∂fF D /∂T 带入）
以减少差分噪声；并对 δT 做敏感性分析。
• 低温区积分分辨率不足：T ≪ Tc 时主要贡献来自 EO(kB T ) 的窄区域。解
决：采用对数网格或分段积分（低能段加密，高能段稀疏），并检查低温幂
律（如 C ∝ T 2 ）是否稳定再用于物理判据。

8

4

热力学求导得比热

本阶段目标是从热力学势或熵的温度导数出发，得到超导态（或有序态）比热曲线 C(T )，
并进一步提取比热跃迁 ∆C、低温幂律等可观测量。常用等价路径包括
C(T ) = T

dS(T )
,
dT

C(T ) = −T

d2 F (T )
,
dT 2

C(T ) =

dU (T )
.
dT

(15)

其中 S 为熵，F 为自由能，U 为内能。若体系存在显著的温度依赖序参量（如能隙
∆(T )）或化学势 µ(T )，严格求导应包含链式法则项，例如


 
 
∂S
d∆
dµ
∂S
∂S
dS
+
=
+
,
(16)
dT
∂T ∆,µ
∂∆ T,µ dT
∂µ T,∆ dT
从而 C(T ) = T dS/dT 会显式包含 d∆/dT 、dµ/dT 等贡献；在弱耦合、自洽求解 ∆(T )
的场景下，上述项往往通过数值上对 S(T ) 直接求导被隐式纳入。

4.1

基于准粒子 DOS 的熵积分与比热求导

文献提取（条目 [6]）给出一条以准粒子态密度（DOS）为核心输入的可复用流程：先
由给定的超导态准粒子 DOS N (E) 与费米–狄拉克分布 fFD (E, T ) 计算熵 Ss (T )，再由
Cs (T ) = T ∂Ss /∂T 得到超导态比热，并以正常态标度 γST 做归一化比较。
定义费米–狄拉克分布
1
,
(17)
fFD (E, T ) = E/(k T )
B
e
+1
其中 E 为准粒子能量，T 为温度，kB 为玻尔兹曼常数。给定每自旋或含自旋简并的
DOS 约定后，可写超导态电子熵（记为 Ss ；下标 s 指 superconducting/state）
Z ∞
h

i
dE N (E) fFD ln fFD + 1 − fFD ln 1 − fFD .
Ss (T ) = −4
(18)
0

式中 N (E) 为超导态准粒子 DOS；因子 4 体现自旋简并与粒子–空穴/支数的计数约定
（应与所用 N (E) 定义保持一致，否则会导致整体系数偏差）。
由熵求比热采用
∂Ss (T )
Cs (T ) = T
.
(19)
∂T
数值实现上，可直接对离散温度网格上的 Ss (T ) 做数值微分得到 Cs (T )。为与实验常用
的线性项系数比较，文献流程采用正常态比热标度 γST 做归一化，报告 Cs /γST 随温度
变化，并读取临界温度 Tc 附近的比热跃迁
∆C
Cs (Tc− ) − Cn (Tc+ )
,
=
C Tc
Cn (Tc+ )

(20)

其中 Cn 为正常态比热（常用近似 Cn ≃ γST T ）。

4.2

低温渐近与节点信息

当能隙具有点节点并导致低能 DOS 呈幂律行为时，低温比热呈现幂律标度。文献流程
指出：点节点各向异性能隙导致低能 DOS 的特定行为，从而得到低温
Cs (T ) ∝ T 2

(T ≪ Tc ),

(21)

该结论可作为数值计算曲线在低温端的自洽检查与物理判据（与线节点常见的 C ∝ T 2
或 T 标度区分需结合具体维度与散射模型；实际应以所用 N (E) 的低能展开为准）。
9

4.3

与自洽能隙求解的接口

在以 N (E) 为输入的框架中，N (E) 常由角向形式因子 f (k̂)、能隙 ∆(T ) 以及散射/能
带参数共同决定。文献提取（[6]）说明其通过弱耦合自洽方程求解 ∆(T ) 以支撑 DOS
的计算，并在自洽方程中施加能量/频率截断 E0 。因此，本阶段与前一阶段（求解 ∆(T )
与构造 N (E)）的关键接口在于：一旦 N (E; T ) 或在每个 T 上由 ∆(T ) 构造的 N (E) 给
定，即可通过式 (18)–(24) 生成 C(T )。
推荐工具
常用工具：
R∞
• 数值积分（ 0 dE）：用于计算熵 Ss (T )（见式 (18)），建议采用自适应
Gauss–Kronrod 或高精度 Simpson，并显式检查能量上限截断的收敛性。
• 数值求解器（自洽方程）：用于求解 ∆(T ) 并生成 N (E)；可用牛顿法/割线
法/拟牛顿法，并在迭代中施加能量/频率截断 E0（文献 [6] 报告使用截断但
未给出软件细节）。
• 数值微分：用于由离散的 Ss (T ) 得到 Cs (T ) = T dSs /dT ；建议采用中心差
分并配合平滑或样条拟合以抑制数值噪声放大。
关键参数
关键参数与调优建议：
• 能 量 积 分 上 限 Emax ： 实 际 计 算 以 E ∈ [0, Emax ] 代 替 [0, ∞)。 推 荐 取
Emax (20–50) kB Tmax ；若 N (E) 在高能端存在结构或引入截断 E0 ，需额
外做 Emax 递增收敛测试。
• 自洽截断 E0 ：弱耦合自洽求解 ∆(T ) 时使用（[6]）。建议将 E0 与配对相互
作用的能标一致，并测试 Tc 与 ∆(0) 对 E0 的敏感性以避免截断主导的伪效
应。
• 温度网格 {Ti }：
为稳定计算 dS/dT ，
推荐在 Tc 附近加密
（例如 T /Tc ∈ [0.8, 1.0]
取更小步长），低温端亦需足够密以分辨幂律（例如检查 C/T 2 是否趋于常
数）。
• 平滑与噪声控制：若 N (E) 来自数值求解或插值，建议对 N (E) 或 S(T ) 做
温和正则化（如样条/局部回归）
，并用不同平滑强度对 ∆C 与低温幂律的影
响做不确定度评估；文献 [6] 未报告对数值噪声与积分误差的专门质控步骤，
复现实作中应补上。

10

案例研究
案例（基于 [6] 的流程复现要点）：
• 输入：两种相（A/B）各自的角向形式因子 f 与相应计算得到的准粒子 DOS
N (E)，并给定正常态比热标度 γST 。
• 计算：用式 (18) 得到 Ss (T )，再用式 (24) 得到 Cs (T )，输出归一化曲线
Cs /γST 。
• 读出物理量：在 Tc 处读取比热跃迁 ∆C/C；在低温区验证点节点情形给出
的幂律 Cs ∝ T 2 。

5

极限/数值评估与特征量

本阶段目标是在给定序参量结构（例如点节点各向异性能隙）或由自洽方程得到的能隙
∆(T ) 支撑下，基于准粒子态密度（density of states, DOS）N (E) 对热力学量进行数值
评估与极限展开：提取低温激活律/幂律、临界温度处比热跳变 ∆C、以及（可扩展地）
磁场/杂质标度关系。以下主要方法来自于以 DOS 为输入的热力学积分框架（见提取
文献 86775778）。

5.1

以 DOS 为输入的熵–比热数值框架

设超导态的准粒子 DOS 为 N (E)（单位：每能量每体积或每自旋通道，具体归一化由
模型约定），费米–狄拉克分布为
fFD (E, T ) =

1
eE/(kB T ) + 1

,

(22)

其中 E ≥ 0 为准粒子能量，kB 为玻尔兹曼常数，T 为温度。以文献中的约定（包含自
旋与粒子/空穴简并的整体因子），超导态熵写为
Z ∞
h

i
Ss (T ) = −4
dE N (E) fFD ln fFD + 1 − fFD ln 1 − fFD ,
(23)
0

其中 Ss 为单位体积熵（或按文中归一化的熵密度）。比热由热力学关系
Cs (T ) = T

∂Ss (T )
∂T

(24)

给出。为便于与实验或正常态线性比热系数比较，引入正常态标度参数 γST （文献中作
为正常态比热系数的标度），并计算归一化曲线
Cs (T )
γST T

或

Cs (T )
.
γST

(25)

在 T = Tc 处定义比热跃迁（跳变）
∆C
∆C
∆C
≡
或
,
+
C
Cn (Tc )
γST Tc

∆C ≡ Cs (Tc− ) − Cn (Tc+ ),

(26)

其中 Cn 为正常态比热（常取 Cn ≃ γST T 作为近 Tc 的主导项；具体取法取决于材料背
景项处理）。
11

5.2

低温极限：点节点导致的幂律

当能隙在费米面上具有点节点（point nodes）时，低能 DOS 常呈幂律
N (E) ∝ E 2 ,

(E ≪ ∆0 ),

(27)

其中 ∆0 为零温能隙量级。将式 (27) 代入式 (23)–(24) 的低温近似可得超导态比热的幂
律标度
Cs (T ) ∝ T 2 ,
(T ≪ Tc ),
(28)
这为区分点节点（T 2 ）与线节点（常见 T ）或全能隙激活律（∼ e−∆/T ）提供了直接的
特征量提取途径。

5.3

能隙自洽求解对 DOS 的支撑与数值实现要点

在弱耦合自洽框架中，可通过能隙方程得到 ∆(T ) 并据此构造 N (E)，再进入熵/比热积
分。文献实现强调：
• 输入：计算得到的 N (E)，以及不同相（例如 A/B 相）的角向形式因子 f (k̂)（用
于描述各向异性能隙结构）。
• 积分：能量积分取 E ∈ [0, ∞)；在自洽方程中施加能量/频率截断 E0（弱耦合常见
做法，用于限制相互作用有效带宽）。
• 方法：数值积分得到 Ss (T )、Cs (T ) 曲线，并从 Tc 邻域读取 ∆C；同时在低温端
检验是否满足式 (28) 的幂律。
上述流程的关键在于（i）N (E) 的低能行为是否与给定节点结构一致；
（ii）T 导数（式
(24)）的数值稳定性；（iii）截断 E0 与能量网格对 ∆(T ) 与热力学量的收敛性影响。
推荐工具
常用工具：
• 数值积分与数值求解器（文献 86775778 报告使用）
：用于计算 ∆(T )、Ss (T )
与 Cs (T ) 曲线并提取 ∆C；软件/平台未在文中注明（实现上可用自编积分
器或通用数值库）。
• 自洽迭代模块：用于弱耦合能隙方程求解（输入截断 E0 与形式因子 f (k̂)，
输出 ∆(T ) 或用于构造 N (E) 的参数）。
• 温度扫描与插值工具：用于在 Tc 附近细化温度步长并稳定读取 Cs (Tc− )、估
计 Cn (Tc+ )。

12

关键参数
关键参数推荐与调优建议：
• 能量截断 E0 ：作为弱耦合自洽的高能截断，应取满足 E0 ≫ ∆0 且结果对
E0 变化不敏感的区间；建议做 E0 扫描并以 ∆(T )、∆C 的收敛为准。
• 能量网格 {Ei }：低能端需加密以解析 N (E) ∝ E 2 的幂律与低温 Cs ∝ T 2 ；
高能端需覆盖至若干倍 kB Tc （或至 E0 ）以保证熵积分尾部收敛。
• 温度步长 ∆T ：在 T ≪ Tc 用对数步长捕捉幂律区间；在 T ≈ Tc 使用更细
线性步长以减少数值微分误差并稳定读取 ∆C。
• 归一化标度 γST ：用于比较不同相或不同参数组的曲线，建议与正常态 Cn /T
的实验或模型值一致；若包含背景项，需明确是否将其并入 γST 或从 Cn 中
扣除。
案例研究
案例（文献 86775778 的典型输出与可读特征量）：
• 以给定的准粒子 DOS N (E) 与相应的角向形式因子 f (k̂)（A/B 相）为输入，
数值积分得到 Ss (T )，再由 Cs = T (∂Ss /∂T ) 生成 Cs (T ) 全曲线。
• 在 Tc 处读取比热跃迁 ∆C 并给出归一化量（如 ∆C/(γST Tc ) 或 ∆C/Cn ），
用于区分不同相的热力学可观测差异。
• 在低温端利用点节点诱导的低能 DOS 行为推得 Cs (T ) ∝ T 2 ，作为节点类型
的幂律判据（与全能隙激活律相对照）。

6

结论

本文围绕超导态比热的计算与表征，总结了不同相结构、节点类型以及磁场/杂质等低
能尺度对热力学量的控制作用。主要结论如下：
1. 比热跃迁与相结构密切相关：以正常态标度 C = γS T 表示时，A 相在 Tc 处给出
∆C/C ≃ 0.93，B 相在 Tc 处给出 ∆C/C ≃ 1.20，显示不同配对对称性/能隙各向
异性可显著改变跃迁幅度。
2. 低温渐近由节点类型主导：对于具有点节点的各向异性两相，低温比热均呈幂律
行为 Cs ∝ T 2 ，体现了低能准粒子态密度的非指数型增长。
3. 在弱耦合近似下，临界温度满足
2γ
E0 e−1/λ = 1.136 E0 e−1/λ ,
Tc =
π

(29)

该关系为将微观耦合常数 λ 与宏观热力学跃迁温度联系起来提供了基准。
4. 允许多分量序参量时可能出现两个相变温度：Tc（∆1 开始非零）与 Tc1（∆2 开始
非零）
，从而在 C(T ) 中对应出现两次跳变；当 ∆1 与 ∆2 同时非零时，dx2 −y2 + idxy
混合态可变为无节点（全能隙），从而改变低温热力学由幂律向指数抑制的趋势
（尽管本工作流未给出具体指数系数）。
13

5. 在考虑杂质展宽尺度 γ0 与磁场诱导能标 EH 时，比热的温区行为可被分段刻画：
在 T ≪ max(γ0 , EH ) ≪ ∆0 的低温极限，
C(T, H) ≃

π2
N (0; H) T,
3

γel =

π2
N (0; H),
3

(30)

呈现类金属的线性温度项；而在 γ0 , EH ≪ T ≪ ∆0 的中温区间，
C(T, H) ≃ N0

9ζ(3) 2
T ,
2

(31)

反映了节点激发主导下的幂律标度。
上述结果仍存在若干局限性：其一，部分关键数值（例如两次跃迁对应的具体跳变
量、全能隙混合态的指数衰减系数等）在现有工作流摘录中缺失，限制了与实验的定量
对照；其二，多数推导依赖弱耦合、准经典或低能有效描述，并采用明确的能标分离
（如 T ≪ ∆0 、γ0 , EH ≪ T 等），在强耦合、复杂多带、显著各向异性费米面或临界涨落
增强的体系中适用性需进一步检验；其三，关于计算流程本身仍缺乏统一的标准化描述
（输入输出规范、误差传播、参数拟合策略与可复现实验对接方式），不利于跨材料、跨
方法的系统比较。
未来工作可沿以下方向推进：第一，推动方法标准化与可复现流程建设，将比热计
算中的假设层级（能隙模型、散射模型、磁场处理、数值积分精度等）模块化，并建立
统一的基准测试（benchmark）以对齐不同实现；第二，增强跨领域适用性，将单带近
似推广至多带与自旋轨道耦合显著体系，并纳入强耦合修正、各向异性电子-声子/自旋
涨落介导机制，以及与输运、穿透深度、NMR 等多实验量的联合约束；第三，提升自
动化能力，发展从实验 C(T, H) 数据到候选配对对称性与参数区间的自动反演流程（包
含不确定度量化与模型选择），从而更高效地区分节点结构、识别多相变并实现材料筛
选与机制判别的闭环优化。

References
[1] Los Angeles,CA90089-0484, USA. Impurity Scattering in f-wave Superconductor
UPt3 . Unknown, Unknown.
[2] V. G. Kogan, C. Martin, † and R. Prozorov ‡ . Superfluid density and specific heat
within self-consistent scheme for two-band superconductor. Unknown, Unknown.
[3] Partha Goswami, Manju Rani, and Avinashi Kapoor. Fermi pockets and quantum
oscillations in specific heat of YBCO in the presence of disorder. Unknown, Unknown.
[4] Qingshan Yuan, Hong-Yi Chen, H. Won, S. Lee, K. Maki, P. Thalmeier, and C. S.
Ting. Impurity effects on s+g-wave superconductivity in borocarbides Y(Lu)Ni�B�C.
Unknown, Unknown.
[5] Tatsuki Hashimoto, Keiji Yada, Ai Yamakage, Masatoshi Sato and Yukio Tanaka.
Bulk Electronic State of Superconducting Topological Insulator. Unknown, Unknown.
[6] David Parker, Kazumi Maki, and Stephan Haas. Anisotropic superconductivity in
PrOs4 Sb12 . Unknown, Unknown.
[7] Giuseppe G.N. Angilella, Renato Pucci, Fabio Siringo. Sharp k-space features in the
order parameter within the interlayer pair-tunneling mechanism of high- Tc superconductivity. Unknown, Unknown.
14

[8] Departamento de Fisica, UFAM, Av. Rodrigo Otavio 3000, Japiim, 69077-000 Manaus, AM, Brazil. Impurity effects on the d-density wave and superconductivity phase
of cuprates. Unknown, 2012.
[9] School of Physics, Georgia Institute of Technology, Atlanta Georgia 30332. F-wave
versus P-wave Superconductivity in Organic Conductors. Unknown, Unknown.
[10] Qingshan Yuan, 1 Xin-Zhong Yan, 1,2 and C. S. Ting 1 . s-Wave-Like excitation in
the superconducting state of electron-doped cuprates with d-wave pairing. Unknown,
Unknown.
[11] Department of Physics, University of Florida, Gainesville, FL 32611, USA.. Vortex Contribution to Specific Heat of Dirty d-Wave Superconductors: Breakdown of
Scaling. Unknown, Unknown.
[12] V. Jovanović, R. Zikic , L. Dobrosavljević-Grujić. Pairing in planar organic superconductors. Physica C, 2005.
[13] Here, N (0) is the total density of states at the Fermi level per one spin. The order
parameter was taken in the form ∆(r, T ; kF ) = Ψ(r, T )Ω(kF ) where Ω(kF ) describes
the variation of ∆ along the Fermi surface and is normalized to have the average over
the whole surface ⟨Ω2 ⟩ = 1 .. Scaling relations in anisotropic superconductors with
strong pair-breaking. Unknown, Unknown.
[14] YU Ya-bin, ZHANG Li-yuan and CHEN Chang-feng. EFFECTIVE MASS AND
SPECIFIC HEAT IN THE PERIODIC ANDERSON MODEL WITH NEGATIVE
U CENTERS. Physica C, 1988.
[15] B. Uchoa 1 , G. G. Cabrera 1 , and A. H. Castro Neto 2 . Nodal liquid and s-wave
superconductivity in transition metal dichalcogenides. Unknown, Unknown.
[16] Instituto de Física Teórica, Universidade Estadual Paulista,. Phase transition from
a dx2 −y2 to dx2 −y2 + idxy superconductor*. Unknown, Unknown.
[17] Department of Physics and Astronomy, University of Southern California, Los Angeles, CA. Anisotropic s-wave superconductivity in MgB2 .Unknown, U nknown.

15

