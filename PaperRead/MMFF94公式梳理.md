## MMFF94_公式梳理

### 目标函数

#### 总能量式子

$$
\begin{align}
E_{\mathrm{MMFF}} =& \sum EB_{ij} + \sum EA_{ijk} + \sum EBA_{ijk} + \sum \mathrm{EOOP}_{ijk;l} \\
&+ \sum E\mathrm{T}_{ijkl} + \sum \mathrm{EvdW}_{ij} + \sum EQ_{ij}
\end{align}
\tag{1}
$$

#### 键伸缩项 （Bond Stretching, Paper III）

$$
EB_{ij} = 143.9325 \frac{k_{b,IJ}}{2} \Delta r_{ij}^{2}
\times (1 + c_s \Delta r_{ij} + 7/12c_s^{2} \Delta r_{ij}^{2})
\\
\Delta r_{ij}=r_{ij}-r_{IJ}^{0}\\
c_{s}=-2 \text{\AA}^{-1}
\tag{2}
$$

| Param | $k_{b,IJ}$        | $r_{ij}$     | $c_s$       |
| ----- | ----------------- | ------------ | ----------- |
| Unit  | $\text{md / \AA}$ | $\text{\AA}$ | $\text \AA$ |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 6_MMFFBOND.csv | III   | III   | 578  |

| BT   | MTYPE_I | MTYPE_J | KB    | R0    | SOURCE |
| ---- | ------- | ------- | ----- | ----- | ------ |
| 0    | 1       | 1       | 4.258 | 1.508 | C94    |
| 0    | 1       | 2       | 4.539 | 1.482 | C94    |

> [!Note]
>
> 遍历方法：按照邻接扩展表遍历，对第i个原子只遍历编号比i大的相邻节点。有几个氢原子直接乘以几。

> [!Warning]
>
> 未知参数$r_{ij}$，可能为优化目标。



#### 角弯曲项 (Angle Bending, Paper III)

1、多数情况下：
$$
EA_{ijk} = 0.043844 \frac{k_{a,IJK}}{2} \Delta \vartheta_{ijk}^{2} (1 + c_b \Delta \vartheta_{ijk})\tag{3}\\
 \Delta \vartheta_{ijk} = \vartheta_{ijk} - \vartheta_{IJK}^0 \\
c_{b}=-0.007 \text{deg}^{-1}
$$
> [!Caution]
>
> 注意deg转rad

2、线性或近线性键角，例如炔烃、腈、异腈、叠氮化物和重氮化合物中出现的键角：
$$
EA_{ijk} = 143.9325 \, k a_{ijk} * (1 + \cos \vartheta_{ijk}) \tag{4}
$$

| Param | $k_{a,IJK}$             | $\vartheta$  | $c_b$             |
| ----- | ----------------------- | ------------ | ----------------- |
| Unit  | $\text{md \AA /rad}^2 $ | $\text{deg}$ | $\text{deg}^{-1}$ |

| 文件          | 论文 | 表格 | 页码 |
| ------------- | ---- | ---- | ---- |
| 8_MMFFANG.csv | III  | IV   | 579  |

| ANGLE_TYPE | MTYPE_I | MTYPE_J | MTYPE_K | KA_IJK | THETA_0 | COMMENT_ORIGIN     |
| ---------- | ------- | ------- | ------- | ------ | ------- | ------------------ |
| 0          | 0       | 1       | 0       | 0      | 108.9   | 0:*-1-* MMFF94 DEF |
| 0          | 1       | 1       | 1       | 0.851  | 109.608 | C94                |
| 0          | 1       | 1       | 2       | 0.736  | 109.445 | C94                |
| 0          | 1       | 1       | 3       | 0.777  | 107.517 | C94                |

> [!Note]
>
> 遍历方法：按照邻接扩展表顺序遍历，无需标记数组。中心原子i，设相邻原子p个，则有$C_m^2=\frac{m(m-1)}{2}$个EA。可与EBA一起计算。
>
> 

> [!Warning]
>
> 未知参数$\vartheta_{ijk}$，可能为优化目标



#### 伸缩-弯曲相互作用  (Stretch-Bend Interactions, Paper III)

$$
\mathrm{EBA}_{ijk} = 2.51210(k_{ba,IJK} \, \Delta r_{ij} + k_{ba,KJI} \, \Delta r_{kj}) \, \Delta \vartheta_{ijk}
\tag{5}
$$

当$k$为正数时，$\vartheta$与$r_{ij}$呈反比。

| Param | $k_{ba,IJK}$      | $r_{ij}$     | $\vartheta$       |
| ----- | ----------------- | ------------ | ----------------- |
| Unit  | $\text{md / rad}$ | $\text{\AA}$ | $\text{deg}^{-1}$ |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 9_MMFFSTBN.csv | III   | V     | 581  |

| SBT  | MTYPE_I | MTYPE_J | MTYPE_K | KBA_IJK | KBA_KJI | SOURCE |
| ---- | ------- | ------- | ------- | ------- | ------- | ------ |
| 0    | 1       | 1       | 1       | 0.206   | 0.206   | C94    |
| 0    | 1       | 1       | 2       | 0.136   | 0.197   | C94    |
| 0    | 1       | 1       | 3       | 0.211   | 0.092   | C94    |
| 0    | 1       | 1       | 5       | 0.227   | 0.07    | C94    |

> [!Note]
>
> 遍历方法同上。

> [!Warning]
>
> 未知参数$r_{ij}$和$\vartheta$，可能为优化目标

> [!Caution]
>
> vartheta注意deg转为rad

#### 三配位中心的面外弯曲 (Out-of-Plane Bending at Tricoordinate Centers, Paper III)

$$
\mathrm{EOOP}_{ijk;l} = 0.043844 \frac{k_{oop,IJK:L}}{2} \chi_{ijk;l}^{2}
\tag{6}
$$

| Param | koop         | $\chi$ |
| ----- | ------------ | ------ |
| Unit  | md Å/rad$^2$ | deg    |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 11_MMFFOOP.csv | III   | VI    | 582  |

| MTYPE_I | MTYPE_J | MTYPE_K | MTYPE_L | KOOP  | SOURCE          |
| ------- | ------- | ------- | ------- | ----- | --------------- |
| 0       | 2       | 0       | 0       | 0.02  | *-2-*-* C94 DEF |
| 1       | 2       | 1       | 2       | 0.03  | C94             |
| 1       | 2       | 2       | 2       | 0.027 | C94             |
| 1       | 2       | 2       | 3       | 0.026 | C94             |

> [!Note]
>
> 遍历方法：当且仅当原子i的AllBondCount=3时才计算该项。

> [!Warning]
>
> 未知参数$\chi$

#### 扭转相互作用项 (Torsion Interactions, Paper IV)

$$
ET_{ijkl} = 0.5\,\bigl[V_{1}(1+\cos\Phi) + V_{2}(1-\cos 2\Phi) + V_{3}(1+\cos 3\Phi)\bigr]\tag{7}
$$

| Param | $\Phi$ | V        |
| ----- | ------ | -------- |
| Unit  | deg    | kcal/mol |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 12_MMFFTOR.csv | IV    | V     | 602  |

| TT   | ATOM_TYPE_I | ATOM_TYPE_J | ATOM_TYPE_K | ATOM_TYPE_L | V1    | V2     | V3    | SOURCE————————    |
| ---- | ----------- | ----------- | ----------- | ----------- | ----- | ------ | ----- | ----------------- |
| 0    | 0           | 1           | 1           | 0           | 0     | 0      | 0.3   | C94 0:*-1-1-* Def |
| 5    | 0           | 1           | 1           | 0           | 0.2   | -0.8   | 1.5   | C94 5:*-1-1-* Def |
| 0    | 1           | 1           | 1           | 1           | 0.103 | 0.681  | 0.332 | C94               |
| 5    | 1           | 1           | 1           | 1           | 0.144 | -0.547 | 1.126 | C94               |

> [!Warning]
>
> 未知参数$\Phi$

> [!Note]
>
> 遍历方法：使用邻接扩展表，增加大小为$1\times n$的标记数组rec（初值为false，表示为遍历）。对于NH表按Seq遍历。设当前为原子i，设相邻原子为s，如果rec[s.seq]=false，则计算is扭转。否则跳过。计算时，若i 除s外有p个未遍历相邻原子，s除i外有q个为遍历相邻原子，则计算次数为$p\times q$。
>
> 遍历方法:

#### 范德华相互作用 (Van Der Waals Interactions, Paper II)

当原子 i 与 j 属于不同分子或被三根及以上化学键隔开时，才计算范德华相互作用 。
$$
EvdW = \varepsilon_{IJ}\, \left(\frac{1.07R_{IJ}^{*}}{R_{ij}+0.07R_{IJ}^{*}}\right)^{7}\, \left(\frac{1.12{R_{IJ}^{*}}^7}{R_{ij}^{7}+0.12{R_{IJ}^{*}}^7} - 2\right)
\tag{8}
$$
参数解释：
$$
\begin{gather}
R_{II}^{*} = A_{I} \alpha_{I}^{1/4}\tag{9}\\
R_{IJ}^{*} = 0.5(R_{II}^{*} + R_{JJ}^{*})(1 + B(1 - \exp(-\beta \gamma_{IJ}^{2})))\tag{10}\\
\gamma_{IJ} = (R_{II}^{*} - R_{JJ}^{*})/(R_{II}^{*} + R_{JJ}^{*})\tag{11}\\
\epsilon_{IJ} = \frac{181.16 G_{I} G_{J} \alpha_{I} \alpha_{J}}{(\alpha_{I}/N_{I})^{1/2} + (\alpha_{J}/N_{J})^{1/2}} \frac{1}{R_{IJ}^{*6}}\tag{12}
\end{gather}
$$
方程修正：

1、一般$B = 0.2$，$\beta = 12$

2、若原子$I$或$J$为极性氢，取$B=0$（暂不考虑）

3、若该相互作用被归类为给体–受体氢键，则将$R_{IJ}^{*} \times DARAD$，将$\varepsilon_{ij} \times DAEPS$。其中$DARAD=0.8$，$DAEPS=0.5$。（暂不考虑）

当前理解：

1、将键中的两个原子当作是一个基团，即当原子i属于MMFF类型I（Symbol）时，则可以计算$R_{II}$

2、原子i与原子j各自从属于不同MMFF类型，且间隔三根键及以上。

| Param | R            | $\varepsilon$ | $R^*$        | N    | G    | $\alpha$ | A    |
| ----- | ------------ | ------------- | ------------ | ---- | ---- | -------- | ---- |
| Unit  | $\text{\AA}$ | kcal/mol      | $\text{\AA}$ |      |      |          |      |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 13_MMFFVDW.csv | II    | III   | 529  |

| MTYPE | ALPHA_I | N_I  | A_I  | G_I   | DA   | SYMBOL | ORIGIN |
| ----- | ------- | ---- | ---- | ----- | ---- | ------ | ------ |
| 1     | 1.05    | 2.49 | 3.89 | 1.282 | -    | CR     | E94    |
| 2     | 1.35    | 2.49 | 3.89 | 1.282 | -    | C=C    | E94    |
| 3     | 1.1     | 2.49 | 3.89 | 1.282 | -    | C=O    | E94    |
| 4     | 1.3     | 2.49 | 3.89 | 1.282 | -    | CSP    | E94    |

> [!Note]
>
> 遍历方法：使用邻接扩展表。设分子共有n个原子，当前为第$i$个原子，设置一个$1\times n$的标记数组rec（初值为False）。以$i$为根节点，进行广度优先搜索，设置搜索深度$\leq$2。将所有遍历到的节点的rec值变为True。然后令j=i...n-1开始循环，若rec[j]=False，则计算$E(i,j)$并加和。
> $$
> \sum_i\sum_jE(i,j)(i\leq j)
> $$
> 

> [!Warning]
>
> 该方程组未知数为$R_{ij}$，理论上可由键角$\vartheta $和相邻原子距离$r_{ij}$导出。

#### 静电相互作用项 (Electrostatic Interactions, Paper II)

要求间隔三根及以上化学键，同范德华作用。
$$
EQ_{ij} =332.0716 \frac{q_i q_j }{ D (R_{ij} + \delta)^n}\tag{13}
$$

$$
q_{i} = q_{I}^{0} + \sum_{K} \omega_{KI}\tag{14}\\
\omega_{IK}=-\omega_{KI}\\
\delta=0.05\text{\AA}
$$

通常取介电常数$D=1$指数 $n=1$，但也支持使用随距离变化的介电常数（$n=2$）

k为i的相邻原子，先归类为MMFFType，然后再通过下表给出求和结果。

| Param | R            |      |      |
| ----- | ------------ | ---- | ---- |
| Unit  | $\text{\AA}$ |      |      |

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 14_MMFFCHG.csv | II    | IV    | 530  |

| CLASS | MTYPE_1 | MTYPE_2 | BCI ($\omega$) | SOURCE |
| ----- | ------- | ------- | -------------- | ------ |
| 0     | 1       | 1       | 0              | #C94   |
| 0     | 1       | 2       | -0.1382        | C94    |
| 0     | 1       | 3       | -0.061         | #C94   |
| 0     | 1       | 4       | -0.2           | #X94   |

==未知参数$R_{ij}$，可以通过$\vartheta$和$r_{ij}$导出==

### 其他补充资料

#### MMFF94_TYPE

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 1_MMFFSYMB.csv | I     | III   | 506  |

| SYMBOL | MTYPE | DEFINITION              |
| ------ | ----- | ----------------------- |
| CR     | 1     | ALKYL CARBON/ SP3       |
| C=C    | 2     | VINYLIC CARBON/ SP2     |
| CSP2   | 2     | GENERIC SP2 CARBON      |
| C=O    | 3     | GENERAL CARBONYL CARBON |
| C=N    | 3     | SP2 CARBON IN C=N       |

#### ATOM_TYPE

当为小写ijkl时直接使用原子序数。当为大写时同MMFF_TYPE

#### H_MTYPE

| File           | Paper | Table | Page |
| -------------- | ----- | ----- | ---- |
| 1_MMFFSYMB.csv | II    | III   | 529  |

| PARANT_MSYMBOL | HYDROGEN_MSYMBOL |
| -------------- | ---------------- |
| CR             | HC               |
| C=C            | HC               |
| CSP2           | HC               |
| C=O            | HC               |
| C=N            | HC               |

#### 初始三维坐标参考值



### 优化方法

**拟牛顿法 (Quasi-Newton Methods)**

拟牛顿法不直接计算海森矩阵，而是通过迭代构造一个近似的海森矩阵。常见的算法有BFGS（Broyden–Fletcher–Goldfarb–Shanno），它在计算复杂度和收敛性上有良好的折衷。

- **优点**：比牛顿法计算量少，收敛速度较快。
- **缺点**：可能需要更多的内存来存储近似的海森矩阵。

**MATLAB实现**：

MATLAB自带的 `fminunc` 函数可以使用拟牛顿法来解决无约束优化问题，通常选用BFGS方法。

```matlab
options = optimset('GradObj', 'on', 'MaxIter', 1000);
[x, fval] = fminunc(@(x) f(x), initial_guess, options);
```

### MMFF结果验证

见Paper II. Table V（P533）。其中给出了一些常见分子的MMFF计算结果，可用于验证。



### 其余补充修正

见Paper V 





### 问题

ATYPE和MTYPE是不是一个TYPE？

$q_I^{0}$怎么确定？单位是什么？

---

### 补充参考资料

统一变量名之前的代码如下，不可删除！！！

```python
# 格式：("文件标题","表格列名",扩展列数目) 扩展列数目为-1则列数未知,0表示无扩展列
ori_title_list = [
    ("1. MMFFSYMB.PAR", "SYMBOL,TYPE,DEFINITION+"),
    ("2. MMFFDEF.PAR", "SYMBOL,TYPE_1,TYPE_2,TYPE_3,TYPE_4,TYPE_5,DEFINITION"),
    ("3. MMFFAROM.PAR", "OLD_TYPE,AROM_TYPE,AT_NUM,RING_SIZE,L5,IM_CAT,N5_ANION"),
    ("4. MMFFHDEF.PAR", "PARANT_TYPE,HYDROGEN_TYPE"),
    ("5. MMFFPROP.PAR", "ATOM_TYPE,ASPEC,CRD,VAL,PILP,MLTB,AROM,LIN,SBMB"),
    ("6. MMFFBOND.PAR", "BT,ATOM_TYPE_I,ATOM_TYPE_J,kb,r0,Source"),
    ("7. MMFFBNDK.PAR", "ATOM_i,ATOM_j,r0(ref),kb(ref),Source"),
    ("8. MMFFANG.PAR", "ANGLE_TYPE,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ka_IJK,theta_0,Comment_origin+"),
    ("9. MMFFSTBN.PAR", "SBT,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,kba_IJK,kba_KJI,Source"),
    ("10. MMFFDFSB.PAR", "IR,JR,KR,F_IJK,F_KJI"),
    ("11. MMFFOOP.PAR", "ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ATOM_TYPE_L,koop,Source+"),
    ("12. MMFFTOR.PAR", "TT,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ATOM_TYPE_L,V1,V2,V3,Source+"),
    ("13. MMFFVDW.PAR", "type,alpha_i,N_i,A_i,G_i,DA,SYMBOL,Origin"),
    ("14. MMFFCHG.PAR", "CLASS,TYPE_1,TYPE_2,bci,Source"),
    ("15. MMFFPBCI.PAR", "_,TYPE,pbci,fcadj,Origin/Comment")
]
```

| File | Paper | Table | Page |
| ---- | ----- | ----- | ---- |
|      |       |       |      |
|      |       |       |      |
|      |       |       |      |

```py
ori_title_list = [
    ("1. MMFFSYMB.PAR", "SYMBOL,TYPE,DEFINITION+", "SYMBOL TYPE  DEFINITION"),
    ("2. MMFFDEF.PAR", "SYMBOL,TYPE_1,TYPE_2,TYPE_3,TYPE_4,TYPE_5,DEFINITION",
     "SYMBOL  TYPE   DEFAULT TYPES DEFINITION"),
    ("3. MMFFAROM.PAR", "OLD_TYPE,AROM_TYPE,AT_NUM,RING_SIZE,L5,IM_CAT,N5_ANION",
     "TYPE   TYPE    NUM SIZE  L5   CAT ANION"),
    ("4. MMFFHDEF.PAR", "", "PARANT_TYPE,HYDROGEN_TYPE", "TYPE    TYPE"),
    ("5. MMFFPROP.PAR", "ATOM_TYPE,ASPEC,CRD,VAL,PILP,MLTB,AROM,LIN,SBMB",
     "atype aspec crd val  pilp mltb arom lin sbm"),
    ("6. MMFFBOND.PAR", "BT,ATOM_TYPE_I,ATOM_TYPE_J,kb,r0,Source", "types       kb         r0    Source"),
    ("7. MMFFBNDK.PAR", "ATOM_i,ATOM_j,r0(ref),kb(ref),Source", "species r0(ref) kb(ref)   Source"),
    ("8. MMFFANG.PAR", "ANGLE_TYPE,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ka,theta_0,Comment_origin+",
     "atom types      ka     theta0    Comment/origin"),
    ("9. MMFFSTBN.PAR", "SBT,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,kba_IJK,kba_KJI,Source"),
    "types I, J, K    kbaIJK    kbaKJI   Source",
    ("10. MMFFDFSB.PAR", "IR,JR,KR,F_IJK,F_KJI", " IR   JR   KR   F(I_J,K)  F(K_J,I)"),
    ("11. MMFFOOP.PAR", "ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ATOM_TYPE_L,koop,Source+",
     "MMFF  atom types     koop    Source"),
    ("12. MMFFTOR.PAR", "TT,ATOM_TYPE_I,ATOM_TYPE_J,ATOM_TYPE_K,ATOM_TYPE_L,V1,V2,V3,Source+",
     "atom types       V1      V2      V3     Source"),
    ("13. MMFFVDW.PAR", "type,alpha_i,N_i,A_i,G_i,DA,SYMBOL,Origin",
     "type  alpha-i     N-i       A-i       G-i DA Symb   Origin"),
    ("14. MMFFCHG.PAR", "CLASS,TYPE_1,TYPE_2,bci,Source", "types       bci     Source"),
    ("15. MMFFPBCI.PAR", "_,TYPE,pbci,fcadj,Origin/Comment", " type    pbci      fcadj   Origin/Comment")
]
```

