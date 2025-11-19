### MMFF94计算示例

> [!Note]
>
> 唯一SMILES：CCC

#### 第一步：初始表格信息补足

将CCC的原子做补足，获得原子表格和邻接表。

##### 非氢原子表（NH表）

CANGEN程序中已给出，从0开始编号。

| Seq  | MType | Element | IsAroma | CCount | CHCount | NHBC | CharSym | CharVal |
| ---- | ----- | ------- | ------- | ------ | ------- | ---- | ------- | ------- |
| 0    |       | C       | 0       | 1      | 3       | 1    | 0       | 0       |
| 1    |       | C       | 0       | 2      | 2       | 2    | 0       | 0       |
| 2    |       | C       | 0       | 1      | 3       | 1    | 0       | 0       |

> [!Warning]
>
> 表格基本实现，仅MType的识别方法待定。

##### 氢原子表（H表）

基于上述非氢原子表生成。编号接非氢原子表格。

| Seq  | ParentSeq | MType |
| ---- | --------- | ----- |
| 3    | 0         | 5     |
| 4    | 0         | 5     |
| 5    | 0         | 5     |
| 6    | 1         | 5     |
| 7    | 1         | 5     |
| 8    | 2         | 5     |
| 9    | 2         | 5     |
| 10   | 2         | 5     |

> [!Warning]
>
> 表格未实现，MType识别方法待定。

##### 邻接扩展表（AL表）

AllBondCount=NHCount+CHCount

```
Seq     NHBond           HBonD                  AllBondCount
0       0: --->1         --->3 --->4 --->5      4
1       1: --->0 --->2   --->6 --->7            4
2       2: --->1         --->8 --->9 --->10     4
```

##### 

> [!Warning]
>
> 表格基本实现，仅HBond未实现。

##### 三维坐标表（XYZ表）

表格包含Seq H。初始节点为0。往正向延伸

| Seq  | X            | Y    | Z    |
| ---- | ------------ | ---- | ---- |
| 0    | 0            | 0    | 0    |
| 1    | $r^{0}_{IJ}$ | 0    | 0    |
| 2    |              |      |      |
| ...  |              |      |      |

可以在生成唯一SMILES同时生成。

| File           | Param | Paper | Table | Page |
| -------------- | ----- | ----- | ----- | ---- |
| 6_MMFFBOND.csv | R0    | III   | III   | 578  |

> [!Note]
>
> 生成方法：
>
> 0、按Seq从上往下生成。设共有n个非H原子，引入一个$1\times n$的bool数组rec，初值false，表示是否遍历。引入一个$1\times 3$数组max_xyz。
>
> 1、设当前遍历到第i个非H原子，共有m个原子与原子i相连。
>
> 2、读取m=AllBondCount，确定构型（键角）。确定方法如下：
>
> ​	（1）按照均布原则，使得n个键在空间均布。
>
> ​		e.g. m=2，直线180度；m=3，平面120度；m=4，正四面体109.5度
>
> ​		m=5，直线+垂直三角形平面。
>
> ​	（2）需要获得构象的顺序编号及相对位置。
>
> e.g. m=4，$RX^2+RY^2+RZ^2=0$
>
> 构象表（球）：
>
> | Seq  | R    | theta(deg) | varphi(deg) |
> | ---- | ---- | ---------- | ----------- |
> | 0    | 1    | 0          | 0           |
> | 1    | 1    | 0          | 109.5       |
> | 2    | 1    | 120        | 109.5       |
> | 3    | 1    | 240        | 109.5       |
>
> 构象表（直）：通过公式生成。
>
> | Seq  | x    | y    | z    |
> | ---- | ---- | ---- | ---- |
> |      |      |      |      |
> |      |      |      |      |
>
> 3、遍历数组rec，其中p个原子已标记。优先Seq最小的将原子$s_0$放入直线键位。
>
> 4、以$\overrightarrow{is_0}$方向为Z轴正方向，以单位向量$\overrightarrow{a}=(a_x,a_y,a_z)$为X轴正方向，单位向量$\overrightarrow{b}=(b_x,b_y,b_z)$为Y轴正方向建系。设max_xyz中取最小值方向为$d_0$，取中值方向为$d_1$,最大值方向为$d_2$（$d\in\{x,y,z\}$。可列出如下方程，可唯一解出$\overrightarrow{a}$。
> $$
> \begin{gather}
> \overrightarrow{is_0}\cdot \overrightarrow{a}=0\\
> |\ \overrightarrow{a}\ |=1\\
> a_{d_2}=0\\
> a_{d_0}\times max\_xyz.d_0\geq0
> \end{gather}
> $$
> ​	又可得方程如下，可唯一解出$b$。
> $$
> \begin{gather}
> \overrightarrow{is_0}\cdot \overrightarrow{b}=0\\
> \overrightarrow{a}\cdot \overrightarrow{b}=0\\
> |\ \overrightarrow{b}\ |=1\\
> \ a_{d_1}\times max\_xyz.d_1\geq0\\
> \end{gather}
> $$
> 5、设空间直角坐标系$I$原点$O_I$，空间直角坐标系$II$原点$O_{II}$。点$O_{II}$在系$I$中的坐标为$(o^{I}_x,o^{I}_y,o^{I}_z)$，点$A$在系$I$中的坐标为$(a_{x}^{I},a_{y}^{I},a_{z}^{I})$，系$II$中的坐标为$(a_{x}^{II},a_{y}^{II},a_{z}^{II})$。已知系$II$的单位基向量在系$I$中表示为$x_{II}^I,y_{II}^I,z_{II}^I$。则：
> $$
> \begin{bmatrix} 
> a_x^I \\ a_y^I \\ a_z^I \end{bmatrix} = 
> \begin{bmatrix} x_{II}^I&y_{II}^I&z_{II}^I
> \end{bmatrix} 
> \begin{bmatrix} a_x^{II} \\ a_y^{II} \\ a_z^{II} \end{bmatrix}
> +\begin{bmatrix} o_x^{I} \\ o_y^{I} \\ o_z^{I} \end{bmatrix}
> $$
> 生成新的构象表。
>
> 6、按照最近原则，以此让剩下的原子$s_1,...,s_{p-1}$占用一个键位。
>
> 7、剩下m-s个原子需要确定构型，依据Seq先后放置。注意跟新max_xyz。
>
> 8、将rec[i]标记为true，读取下一个。

综上，可以获得一个初始的三维坐标表。
$$
\vec{N_jN_i}=
$$


##### 补充：环构象预计算

探查环的个数，寻找是否有能够满足所有ref值的解。

#### 第二步：遍历方法确定

一遍主循环完成所有遍历，具体实现见matlab程序。

#### 第三步V1：参数优化

使用xyz做fmincon和粒子群优化均效果不好，原因在于xyz坐标可取值范围过于宽泛，导致效果不佳。

根据分子的伸缩、弯曲、旋转的特性。决定使用初始三维坐标计算方法，加以随机可变量。

##### 键长可变量

键长可变量（$\Delta r_{change,\ ij}$）和设置$ratio_r=\frac{\Delta r_{change,\ ij}}{r_0}$，暂定$ratio \in [1.1,\ 0.9]$。

##### 键偏移角度可变量

$\alpha_{change,ij}$

##### 键正常旋转角度可变量

$\beta_{change,j}$，注意，当且仅当只有一个键确定时有该变量。

#### 第三步V2：参数优化第二版

##### 键正常旋转

1、粗略计算：

​	首先按编号对于每个原子在[0，360]内做计算。设该分子键数为m=4，则每隔a=360/(m-1)/2取一次值，找最小值。

​	这样对于n个分子只需计算6n次。

​	终止条件：某一轮所有的abs(old_alpha-new_alpha)<a

2、精细计算1

​	在alpha+[-a,a]范围内做计算。二分法（黄金分割法）找最小值。

​	终止条件：某一轮所有的abs(new_E-old_E)<1e-4或者abs(new_D-old_D)<1e-2

3、精细计算2（对于m=3，）

​	补充：环丙烷张角60度，sp2杂化键角120度。

​	在beta+[-b,b]范围内计算，黄金分割法找最小值。



##### 键偏移角度

