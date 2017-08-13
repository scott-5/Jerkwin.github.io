---
 layout: post
 title: day2
 categories:
 - 科
 tags:
 - gmx
 math: true
---

# 配置msys2

## 环境变量

`/home/【用户名】/.bashrc`

`export PATH=$PATH:/C/GMX/amber14/bin`
`export AMBERHOME=/C/GMX/amber14`

# PDB文件基础

## 下载pdb文件

类似 [溶菌酶教程](http://jerkwin.github.io/GMX/GMXtut-1/#第一步-准备拓扑)

## 熟悉pdb文件

阅读 [PDB文件格式说明](https://jerkwin.github.io/2015/06/05/PDB%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E8%AF%B4%E6%98%8E)

重点关注

- 链: `CHAIN`
- 缺失原子: `MISSING ATOM`
- 配体信息: `HET`, `HETNAM`
- 二硫键: `SSBOND`

## 处理pdb文件

以`4xp1`为例, 为多巴胺复合物

简单文本处理Notepad2+bash

3. 删除 杂原子: `HETATM`
2. 删除 连接信息: `CONECT`
1. 删除 水分子: `HETATM` `HOH`
2. 抽取 多巴胺: `LDP`

# 多巴胺模拟: 小分子示例

## 获取构型

1. 从pdb中抽取
2. 网络查询
2. 自己构建

确保构型完整. 若pdb中缺失原子, 务必补充完整

## 创建GAFF力场拓扑

### 1. `antechamber`: 计算电荷并生成mol2文件

电荷可使用两种模型:

1. RESP: 基于量化静电势拟合, 计算较慢, 适用于小分子
2. AM1-BCC: 半经验拟合, 快速, 可处理大分子

示例AM1-BCC电荷模型

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>antechamber -i ldp.pdb -fi pdb -o ldp.mol2 -fo mol2 -pf y -c bcc
</pre></div>

如果生成的`ldp.mol2`文件正常, 则说明命令运行成功. 其他文件可删除.

生成的mol2文件, 有时存在不一致之处, 需要改正

<table class="highlighttable"><th colspan="2" style="text-align:left">ldp.mol2</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>@&lt;TRIPOS&gt;MOLECULE
LDP                            <span style="color: #008800; font-style: italic"># 分子名称</span>
   <span style="color: #666666">22</span>    <span style="color: #666666">22</span>     <span style="color: #666666">1</span>     <span style="color: #666666">0</span>     <span style="color: #666666">0</span>  <span style="color: #008800; font-style: italic"># 原子数, 键数, 分子数目</span>
SMALL
bcc

@&lt;TRIPOS&gt;ATOM
 <span style="color: #008800; font-style: italic"># 编号 原子         x         y         z        分子编号 名称     电荷</span>
	  <span style="color: #666666">1</span> C7         -9.1900    2.9910  -26.8180 c3      <span style="color: #666666">708</span> LDP     -0.038100
	  <span style="color: #666666">2</span> C1        -10.3400    3.5770  -26.0270 ca      <span style="color: #666666">708</span> LDP     -0.075300
	  ...
</pre></div>
</td></tr></table>

如果`ATOM`部分中的分子编号和名称与前面的不一致, 需要修改至相同, 否则后续步骤会出错.

我们得到的`ldp.mol2`中两处的分子编号不一致, 所以要将`708`改为`1`, 如下

<table class="highlighttable"><th colspan="2" style="text-align:left">ldp.mol2</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>@&lt;TRIPOS&gt;MOLECULE
LDP                            <span style="color: #008800; font-style: italic"># 分子名称</span>
   <span style="color: #666666">22</span>    <span style="color: #666666">22</span>     <span style="color: #666666">1</span>     <span style="color: #666666">0</span>     <span style="color: #666666">0</span>  <span style="color: #008800; font-style: italic"># 原子数, 键数, 分子数目</span>
SMALL
bcc

@&lt;TRIPOS&gt;ATOM
 <span style="color: #008800; font-style: italic"># 编号 原子         x         y         z        分子编号 名称     电荷</span>
	  <span style="color: #666666">1</span> C7         -9.1900    2.9910  -26.8180 c3        <span style="color: #666666">1</span> LDP     -0.038100
	  <span style="color: #666666">2</span> C1        -10.3400    3.5770  -26.0270 ca        <span style="color: #666666">1</span> LDP     -0.075300
	  ...
</pre></div>
</td></tr></table>

### 2. `parmchk2`: 检查GAFF参数并生成缺失参数文件

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>parmchk2 -i ldp.mol2 -f mol2 -o ldp.frcmod
</pre></div>

### 3. `sleap`: 生成AMBER参数文件和坐标文件

编写`sleap`需要的输入文件, 文件名任意, 如`leap.in`, 内容如下:

<table class="highlighttable"><th colspan="2" style="text-align:left">leap.in</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%">1
2
3
4
5
6
7</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span><span style="color: #AA22FF">source</span> leaprc.ff14SB
<span style="color: #AA22FF">source</span> leaprc.gaff
loadamberparams ldp.frcmod
<span style="color: #B8860B">lig</span><span style="color: #666666">=</span>loadmol2 ldp.mol2
check lig
saveamberparm lig ldp.prmtop ldp.inpcrd
quit
</pre></div>
</td></tr></table>

运行`slaep`命令

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>sleap -f leap.in
</pre></div>

现在`ldp`目录内容如下

### 4. `acpype`: 将AMBER文件转换为GROMACS文件

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>acpype -p ldp.prmtop -x ldp.inpcrd -d
</pre></div>

现在目录内容如下

我们得到了运行GROMACS的坐标文件`ldp_GMX.gro`和拓扑文件`ldp_GMX.top`.

## 真空单分子

为了验证所得的多巴胺分子的拓扑是否正常, 可以运行单个多巴胺分子的NVT模拟. 当然, 如果你经验丰富, 也可以直接模拟更复杂的体系. 但如果出现问题, 排查错误时不易定位.

### 指定盒子大小

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx editconf -f ldp_GMX.gro -d <span style="color: #666666">1</span> -o conf.gro
</pre></div>

### 准备拓扑文件

将`ldp_GMX.top`更名为默认名`topol.top`

### 准备mdp文件

分子很简单, 我们不进行预平衡, 直接做NVT模拟.

<table class="highlighttable"><th colspan="2" style="text-align:left">grompp.mdp</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span><span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 运行控制参数</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
integrator       <span style="color: #666666">=</span> md<span style="color: #008800; font-style: italic">     ; 积分方法, md: 蛙跳; sd: 随机; steep/cg/lbfgs: 能量最小化</span>
dt               <span style="color: #666666">=</span> <span style="color: #666666">1E-3</span><span style="color: #008800; font-style: italic">   ; 积分步长(ps), EM不用</span>
nsteps           <span style="color: #666666">=</span> <span style="color: #666666">10000</span><span style="color: #008800; font-style: italic">  ; 最大积分步数, 默认0表示不限制</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 输出控制参数</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
nstxout                 <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">    ; trr坐标的输出频率(步)</span>
nstvout                 <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">    ; 速度</span>
nstfout                 <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">    ; 力</span>
nstxout<span style="color: #666666">-</span>compressed      <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">    ; xtc压缩坐标的输出频率</span>
compressed<span style="color: #666666">-</span>x<span style="color: #666666">-</span>precision  <span style="color: #666666">=</span> <span style="color: #666666">1000</span><span style="color: #008800; font-style: italic"> ; xtc坐标的精度</span>

nstlog                  <span style="color: #666666">=</span> <span style="color: #666666">1000</span><span style="color: #008800; font-style: italic"> ; 日志文件中能量的输出频率</span>
nstenergy               <span style="color: #666666">=</span> <span style="color: #666666">1000</span><span style="color: #008800; font-style: italic"> ; 能量文件</span>
nstcalcenergy           <span style="color: #666666">=</span> <span style="color: #666666">100</span><span style="color: #008800; font-style: italic">  ; 计算能量的频率, 最好为nstlistt的倍数</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 邻区搜索参数</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
nstlist             <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">      ; 邻区列表更新频率, 0: 用于真空模拟; -1: 自动;</span>
rlist               <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">      ; 邻区列表的截断距离(nm)</span>
cutoff<span style="color: #666666">-</span>scheme       <span style="color: #666666">=</span> Verlet<span style="color: #008800; font-style: italic"> ; 截断方式(Verlet: 粒子截断; group: 电荷组)</span>
ns_type             <span style="color: #666666">=</span> grid<span style="color: #008800; font-style: italic">   ; 邻区搜索方法(grid: 较快; simple: 仅与group联用)</span>

pbc                 <span style="color: #666666">=</span> xyz<span style="color: #008800; font-style: italic">    ; 周期性比较条件: xyz, xy, no: 忽略盒子, 截断与nstlist置零</span>
periodic<span style="color: #666666">-</span>molecules  <span style="color: #666666">=</span> no<span style="color: #008800; font-style: italic">     ; 周期性分子: no, yes</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 静电与范德华</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
rvdw         <span style="color: #666666">=</span> <span style="color: #666666">.9</span><span style="color: #008800; font-style: italic">      ; 范德华截断半径</span>
rcoulomb     <span style="color: #666666">=</span> <span style="color: #666666">.9</span><span style="color: #008800; font-style: italic">      ; 静电截断半径</span>
vdw<span style="color: #666666">-</span>type     <span style="color: #666666">=</span> Cut<span style="color: #666666">-</span>off<span style="color: #008800; font-style: italic"> ; 范德华计算方法</span>
coulombtype  <span style="color: #666666">=</span> PME<span style="color: #008800; font-style: italic">     ; 静电的计算方法</span>

DispCorr     <span style="color: #666666">=</span> No<span style="color: #008800; font-style: italic">   ; 能量/压力的色散长程校正, no: 无; Ener: 能量; EnerPres: 能量和压力</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 温度耦合</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
tcoupl           <span style="color: #666666">=</span> nose<span style="color: #666666">-</span>hoover<span style="color: #008800; font-style: italic"> ; 耦合方法, no: 无; v-rescale: 快速; nose-hoover: 精确</span>
tc<span style="color: #666666">-</span>grps          <span style="color: #666666">=</span> system<span style="color: #008800; font-style: italic">      ; 温度耦合组, 可多个</span>
tau_t            <span style="color: #666666">=</span> <span style="color: #666666">2</span><span style="color: #008800; font-style: italic">           ; 耦合时间常数(ps)</span>
ref_t            <span style="color: #666666">=</span> <span style="color: #666666">298.15</span><span style="color: #008800; font-style: italic">      ; 参考温度(K)</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 压力耦合</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
pcoupl           <span style="color: #666666">=</span> No<span style="color: #008800; font-style: italic">        ; 耦合方法, no: 无, 盒子大小不变; berendsen: 快速; Parrinello-Rahman: 精确</span>
pcoupltype       <span style="color: #666666">=</span> Isotropic<span style="color: #008800; font-style: italic"> ; 耦合类型, isotropic: 各向同性;</span>
<span style="color: #008800; font-style: italic">							 ; semiisotropic: x/y方向各向同性, 与z方向不同, 膜模拟</span>
<span style="color: #008800; font-style: italic">							 ; anisotropic: 各向异性, 盒子可能剧烈变形</span>
<span style="color: #008800; font-style: italic">							 ; surface-tension: 表面张力</span>
tau<span style="color: #666666">-</span>p            <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">         ; 时间常数(ps)</span>
compressibility  <span style="color: #666666">=</span> <span style="color: #666666">4.5E-5</span><span style="color: #008800; font-style: italic">    ; 压缩率(1/bar)</span>
ref<span style="color: #666666">-</span>p            <span style="color: #666666">=</span> <span style="color: #666666">1</span><span style="color: #008800; font-style: italic">         ; 参考压力(bar)</span>

<span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 初始速度</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
gen_vel           <span style="color: #666666">=</span> yes<span style="color: #008800; font-style: italic">    ; no: 使用gro文件的值; yes: 随机</span>
gen_temp          <span style="color: #666666">=</span> <span style="color: #666666">298.15</span><span style="color: #008800; font-style: italic"> ; 随机速度对应的温度</span>
gen<span style="color: #666666">-</span>seed          <span style="color: #666666">=</span> <span style="color: #666666">-1</span><span style="color: #008800; font-style: italic">     ; 随机数种子; -1: 自动确定</span>
</pre></div>
</td></tr></table>

### 执行模拟

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx grompp
gmx mdrun
</pre></div>

查看原子运行是否有明显异常

## 多分子盒子

多巴胺盒子的模拟, 示例如何创建溶剂盒子.

### 创建溶剂盒子

我们使用前面使用的盒子, 将多巴胺分子插入其中

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx insert-molecules -f conf.gro -ci ldp.pdb -o ldp_box.gro -nmol <span style="color: #666666">1000</span> -try 200
</pre></div>

`-nmol`指定要插入的分子数目, 如果要填满, 可设个很大的值
`-try`指定每次插入的尝试次数

在我的机器上运行, 插入了74个多巴胺.

得到的构型可能不够好, 如果需要更合理的构型, 可以试试packmol

### 准备拓扑文件

只须修改多巴胺的数目即可

<table class="highlighttable"><th colspan="2" style="text-align:left">topol.top</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%">1
2
3
4</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span><span style="color: #666666">...</span>
[ molecules ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">; Compound        nmols</span>
 ldp              <span style="color: #666666">75</span>
</pre></div>
</td></tr></table>

### 准备mdp文件

NVT使用前面的即可, NPT则修改压力耦合部分`pcoupl`的值即可

### 执行模拟

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx grompp -c ldp_box.gro -maxwarn 2
</pre></div>

可能遇到警告. 明白意思, 为什么, 修正或忽略

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx mdrun
</pre></div>

运行中可能出错

- 未平衡好, 可先预平衡一下: EM, NVT, NPT
- 时间步长过大, 先改小, 平衡好后再使用大步长
- 刚性水分子无法约束, 先换用柔性水模型

## 单分子水溶液

### 创建水溶液

方法多样:

1. `gmx solvate`
2. `gmx insert-molecules`
3. packmol
4. 其他软件

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx solvate -cp conf.gro -cs -o ldp_wat.gro
</pre></div>

### 准备拓扑文件

我们前面使用的是独立的拓扑文件, 其中没有引用任何其他文件. 对简单体系, 这种文件用起来简单, 但对复杂体系, 容易出问题, 所以我们修改文件这个文件, 将其变为itp文件, 方便引用.

改为itp的方法: 只保留与多巴胺有关的部分

<table class="highlighttable"><th colspan="2" style="text-align:left">ldp.itp</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>[ atomtypes ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb</span>
 c3       c3          <span style="color: #666666">0.00000</span>  <span style="color: #666666">0.00000</span>   A     <span style="color: #666666">3.39967e-01</span>   <span style="color: #666666">4.57730e-01</span><span style="color: #008800; font-style: italic"> ; 1.91  0.1094</span>
 ca       ca          <span style="color: #666666">0.00000</span>  <span style="color: #666666">0.00000</span>   A     <span style="color: #666666">3.39967e-01</span>   <span style="color: #666666">3.59824e-01</span><span style="color: #008800; font-style: italic"> ; 1.91  0.0860</span>
<span style="color: #666666">...</span>

[ moleculetype ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;name            nrexcl</span>
 ldp              <span style="color: #666666">3</span>

[ atoms ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type</span>
     <span style="color: #666666">1</span>   c3     <span style="color: #666666">1</span>   LDP    C7    <span style="color: #666666">1</span>    <span style="color: #666666">-0.038100</span>     <span style="color: #666666">12.01000</span><span style="color: #008800; font-style: italic"> ; qtot -0.038</span>
     <span style="color: #666666">2</span>   ca     <span style="color: #666666">1</span>   LDP    C1    <span style="color: #666666">2</span>    <span style="color: #666666">-0.075300</span>     <span style="color: #666666">12.01000</span><span style="color: #008800; font-style: italic"> ; qtot -0.113</span>
<span style="color: #666666">...</span>

[ bonds ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;   ai     aj funct   r             k</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">2</span>   <span style="color: #666666">1</span>    <span style="color: #666666">1.5130e-01</span>    <span style="color: #666666">2.7070e+05</span><span style="color: #008800; font-style: italic"> ;     C7 - C1</span>
     <span style="color: #666666">1</span>     <span style="color: #666666">10</span>   <span style="color: #666666">1</span>    <span style="color: #666666">1.5350e-01</span>    <span style="color: #666666">2.5363e+05</span><span style="color: #008800; font-style: italic"> ;     C7 - C8</span>
<span style="color: #666666">...</span>

[ pairs ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;   ai     aj    funct</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">6</span>      <span style="color: #666666">1</span><span style="color: #008800; font-style: italic"> ;     C7 - C5</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">7</span>      <span style="color: #666666">1</span><span style="color: #008800; font-style: italic"> ;     C7 - C3</span>
<span style="color: #666666">...</span>

[ angles ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">;   ai     aj     ak    funct   theta         cth</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">2</span>      <span style="color: #666666">4</span>      <span style="color: #666666">1</span>    <span style="color: #666666">1.2063e+02</span>    <span style="color: #666666">5.3421e+02</span><span style="color: #008800; font-style: italic"> ;     C7 - C1     - C2</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">2</span>      <span style="color: #666666">5</span>      <span style="color: #666666">1</span>    <span style="color: #666666">1.2063e+02</span>    <span style="color: #666666">5.3421e+02</span><span style="color: #008800; font-style: italic"> ;     C7 - C1     - C6</span>
<span style="color: #666666">...</span>

[ dihedrals ]<span style="color: #008800; font-style: italic"> ; propers</span>
<span style="color: #008800; font-style: italic">; for gromacs 4.5 or higher, using funct 9</span>
<span style="color: #008800; font-style: italic">;    i      j      k      l   func   phase     kd      pn</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">2</span>      <span style="color: #666666">4</span>      <span style="color: #666666">7</span>      <span style="color: #666666">9</span>   <span style="color: #666666">180.00</span>  <span style="color: #666666">15.16700</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C7-    C1-    C2-    C3</span>
     <span style="color: #666666">1</span>      <span style="color: #666666">2</span>      <span style="color: #666666">4</span>     <span style="color: #666666">14</span>      <span style="color: #666666">9</span>   <span style="color: #666666">180.00</span>  <span style="color: #666666">15.16700</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C7-    C1-    C2-    H2</span>
<span style="color: #666666">...</span>

[ dihedrals ]<span style="color: #008800; font-style: italic"> ; impropers</span>
<span style="color: #008800; font-style: italic">; treated as propers in GROMACS to use correct AMBER analytical function</span>
<span style="color: #008800; font-style: italic">;    i      j      k      l   func   phase     kd      pn</span>
     <span style="color: #666666">2</span>      <span style="color: #666666">6</span>      <span style="color: #666666">5</span>     <span style="color: #666666">15</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C1-    C5-    C6-    H6</span>
     <span style="color: #666666">2</span>      <span style="color: #666666">7</span>      <span style="color: #666666">4</span>     <span style="color: #666666">14</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C1-    C3-    C2-    H2</span>
     <span style="color: #666666">3</span>      <span style="color: #666666">4</span>      <span style="color: #666666">7</span>      <span style="color: #666666">8</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C4-    C2-    C3-    O1</span>
     <span style="color: #666666">3</span>      <span style="color: #666666">5</span>      <span style="color: #666666">6</span>     <span style="color: #666666">16</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C4-    C6-    C5-    H5</span>
     <span style="color: #666666">4</span>      <span style="color: #666666">5</span>      <span style="color: #666666">2</span>      <span style="color: #666666">1</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C2-    C6-    C1-    C7</span>
     <span style="color: #666666">6</span>      <span style="color: #666666">7</span>      <span style="color: #666666">3</span>      <span style="color: #666666">9</span>      <span style="color: #666666">4</span>   <span style="color: #666666">180.00</span>   <span style="color: #666666">4.60240</span>   <span style="color: #666666">2</span><span style="color: #008800; font-style: italic"> ;     C5-    C3-    C4-    O2</span>
</pre></div>
</td></tr></table>

总的拓扑文件

<table class="highlighttable"><th colspan="2" style="text-align:left">topol.top</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span># include <span style="color: #BB4444">&quot;C:\GMX\GMX5.1.4\share\gromacs\top\amber99sb-ildn.ff\forcefield.itp&quot;</span>

# include <span style="color: #BB4444">&quot;ldp.itp&quot;</span>

# include <span style="color: #BB4444">&quot;C:\GMX\GMX5.1.4\share\gromacs\top\amber99sb-ildn.ff\spce.itp&quot;</span>

[ system ]
 ldp

[ molecules ]<span style="color: #008800; font-style: italic"></span>
<span style="color: #008800; font-style: italic">; Compound        nmols</span>
 ldp              <span style="color: #666666">1</span>
 SOL              <span style="color: #666666">550</span>
</pre></div>
</td></tr></table>

__注意__: 由于`ldp.itp`中包含`[ atomtypes ]`部分, 在总拓扑中需要注意各个itp文件的顺序.

### 准备mdp文件

### 执行模拟

## 溶剂化能

分子在水中的溶剂化能, 也就是水合自由能, 可表征分子的水溶性大小.

自由能计算很重要, 却困难, 而且麻烦. 针对不同体系, 有不同的处理方法:

- 小分子使用较准确的方法: `gmx bar`, 须在不同条件下运行多次模拟
- 蛋白质使用较粗略的方法: `MMPBSA`, 须借助外部程序
- 更实用的分析: 相互作用能, 牵引势能

我们以多巴胺溶剂化能为例, 学习自由能的计算

### 创建水盒子并充分平衡

使用前面的水盒子, 作为示例, 没有充分平衡好

### 准备拓扑文件

无须修改

### 准备mdp文件

自由能计算的关键步骤, mdp文件中涉及到的选项

<table class="highlighttable"><th colspan="2" style="text-align:left">grompp.mdp</th><tr><td><div class="linenodiv" style="background-color: #f0f0f0; padding-right: 10px"><pre style="line-height: 125%"> 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39</pre></div></td><td class="code"><div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span><span style="color: #008800; font-style: italic">;=======================================================================</span>
<span style="color: #008800; font-style: italic">; 自由能计算</span>
<span style="color: #008800; font-style: italic">;=======================================================================</span>
free<span style="color: #666666">-</span>energy            <span style="color: #666666">=</span> no
couple<span style="color: #666666">-</span>moltype         <span style="color: #666666">=</span>

sc<span style="color: #666666">-</span>power               <span style="color: #666666">=</span> <span style="color: #666666">1</span>
sc<span style="color: #666666">-</span>alpha               <span style="color: #666666">=</span> <span style="color: #666666">0</span>
sc<span style="color: #666666">-</span>sigma               <span style="color: #666666">=</span> <span style="color: #666666">0.3</span>
sc<span style="color: #666666">-</span>r<span style="color: #666666">-</span>power             <span style="color: #666666">=</span> <span style="color: #666666">6</span>
sc<span style="color: #666666">-</span>coul                <span style="color: #666666">=</span> no

couple<span style="color: #666666">-</span>intramol        <span style="color: #666666">=</span> no

couple<span style="color: #666666">-</span>lambda0         <span style="color: #666666">=</span> vdw<span style="color: #666666">-</span>q
couple<span style="color: #666666">-</span>lambda1         <span style="color: #666666">=</span> vdw<span style="color: #666666">-</span>q
init<span style="color: #666666">-</span>lambda<span style="color: #666666">-</span>state      <span style="color: #666666">=</span> <span style="color: #666666">-1</span>

coul<span style="color: #666666">-</span>lambdas           <span style="color: #666666">=</span>
vdw<span style="color: #666666">-</span>lambdas            <span style="color: #666666">=</span>

fep<span style="color: #666666">-</span>lambdas            <span style="color: #666666">=</span>
mass<span style="color: #666666">-</span>lambdas           <span style="color: #666666">=</span>
bonded<span style="color: #666666">-</span>lambdas         <span style="color: #666666">=</span>
restraint<span style="color: #666666">-</span>lambdas      <span style="color: #666666">=</span>
temperature<span style="color: #666666">-</span>lambdas    <span style="color: #666666">=</span>

init<span style="color: #666666">-</span>lambda            <span style="color: #666666">=</span> <span style="color: #666666">-1</span>
delta<span style="color: #666666">-</span>lambda           <span style="color: #666666">=</span> <span style="color: #666666">0</span>

calc<span style="color: #666666">-</span>lambda<span style="color: #666666">-</span>neighbors  <span style="color: #666666">=</span> <span style="color: #666666">1</span>
init<span style="color: #666666">-</span>lambda<span style="color: #666666">-</span>weights    <span style="color: #666666">=</span>

nstdhdl                <span style="color: #666666">=</span> <span style="color: #666666">50</span>
dhdl<span style="color: #666666">-</span>print<span style="color: #666666">-</span>energy      <span style="color: #666666">=</span> no
separate<span style="color: #666666">-</span>dhdl<span style="color: #666666">-</span>file     <span style="color: #666666">=</span> yes
dhdl<span style="color: #666666">-</span>derivatives       <span style="color: #666666">=</span> yes
dh_hist_size           <span style="color: #666666">=</span> <span style="color: #666666">0</span>
dh_hist_spacing        <span style="color: #666666">=</span> <span style="color: #666666">0.1</span>
</pre></div>
</td></tr></table>

首先, 我们需要指定要计算哪个分子的自由能, `couple-moltype`, 这需要使用索引文件来指定.

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx make_ndx -f ldp_wat.gro
</pre></div>

其次, 使用`integrator = sd`, `pcoupl = Parrinello-Rahman`否则会警告.

<div class="highlight" style="background: #f8f8f8"><pre style="line-height: 125%"><span></span>gmx grompp
</pre></div>

### 执行模拟

理论上, 对每个状态都需要做预平衡, 但对于小分子, 由于比较容易平衡, 所以我们可以使用某一状态的平衡构型(最好是中间状态)直接进行自由能计算即可.

## 伞形采样

计算多巴胺的跨POPC膜的PMF

### 创建POPC双层膜

有多种方法:

1. `gmx genconf`: 堆叠
2. packmol: 自定义性好
3. VMD: 根据已平衡好构型扩展
4. CHARMM-GUI: 功能最强
5. 自编脚本
6. LipidBook下载

所得构型必须进行重新平衡, 性质验证后方可进行成品模拟.

我们使用最简单的VMD

`Extensions`->`Modeling`->`Membrane Builder`

创建30x30的POPC膜, VMD同时会添加水分子

转为gro文件: `gmx  editconf -f 2 2 10`

抽取单个POPC分子, 获取其GAFF拓扑文件, 并整理为itp文件备用

使用外来构型的问题在于, 原子编号顺序与默认的可能不一致

### 准备拓扑文件和mdp文件

关键在于理解伞形采样的实质: 在反应坐标上取一系列相互间有重叠的点, 对每个点进行外加偏离势作用下的模拟, 对所有点的能量数据使用WHAM(加权直方分析方法)计算PMF

模拟流程:

1. 定义反应坐标
2. 反应坐标取点
3. 外加偏离势模拟
4. WHAM分析计算PMF

其中2是难点, 可根据不同体系采用不同方法, 不必照搬, 只要所得构型满足要求即可.

计算多巴胺跨POPC膜的PMF:

1. 反应坐标: 垂直膜方向, 从一端到另一端
2. 取点方法: 牵引, 非平衡, 虚拟键长
3. 偏离势: 简谐势
4. WHAM: `gmx wham`

使用与其他粒子无相互作用, 位置固定的参考原子, 以其为牵引的参考组

抽取指定坐标的点, 以每个点为初始构型, 设牵引速率为零, 做模拟

使用`gmx wham`分析所得数据

# 4xp1模拟: 蛋白示例

## 处理pdb文件

检查构型

使用`spdbv`补充缺失原子

## 创建amber力场拓扑

`gmx pdb2gmx`

## 单独模拟

真空

## 水溶液模拟

加水: `gmx solvate`
离子: `gmx genion`或`gmx insert-molecules`
预平衡: EM, NVT, NPT

# 4xp1-多巴胺: 蛋白配体示例

## 基本流程

1. 准备蛋白
	自行创建/下载蛋白PDB
	PDB数据清洗
	pdb2gmx 获取蛋白拓扑文件
2. 准备配体
	配体构型文件
	配体拓扑文件: GAFF
3. 装配复合物
4. 组合复合物拓扑文件
5. 添加水/离子
6. 预平衡: EM, PR, NVT, NPT
7. 成品模拟 NPT
8. 分析

## 装配复合物

配体结合位置(蛋白配体相对位置)

- PDB中已有结合位置
- 自行决定结合位置: 模拟, 对接, 经验

## 组合复合物拓扑文件

手动编写

## 添加水/离子

- `gmx solvate`
- `gmx insert-molecules`
- `gmx genion`

## 预平衡: EM, PR, NVT, NPT

## 成品模拟: NPT

## 分析

### 平衡判断: RMSD, RMSF, 二级结构

### 相互作用: 结合能

# 蛋白-配体MMPBSA结合能

## 自由能计算方法

- 小分子 准确方法: 自由能微扰(Zwanzig), 接受比例(Bennet, `gmx bar`), 热力学积分(Kirkwood), 非平衡近似(Jarzynski)
- 蛋白 经验方法: MMGBSA, MMPBSA

## MMPBSA

### 理论基础

$$\D G=G_\text{cmp}-G_\text{pro}-G_\text{lig}$$

$$\alg
G &=G_\text{gas}+G_\text{slv} \\
 &= [ E_\text{gas}-TS_\text{gas} ] + [ G_\text{polar}+G_\text{nonpolar} ] \\
 &= [ E_\text{MM}-TS_\text{MM} ] + [ G_\text{PB}+G_\text{surface} ] \\
 &= [ E_\text{MM}-TS_\text{MM} ] + [ G_\text{PB}+\g A+b ] \\
\ealg$$

### 可选程序

皆依赖于GROMACS, APBS

<table id='tab-0'><caption>比较</caption>
<tr>
  <th rowspan="1" colspan="1" style="text-align:center;">程序</th>
  <th rowspan="1" colspan="1" style="text-align:center;">语言</th>
  <th rowspan="1" colspan="1" style="text-align:center;">编译</th>
  <th rowspan="1" colspan="1" style="text-align:center;">系统</th>
  <th rowspan="1" colspan="1" style="text-align:center;">效率</th>
  <th rowspan="1" colspan="1" style="text-align:center;">能量分解</th>
  <th rowspan="1" colspan="1" style="text-align:center;">丙氨酸突变</th>
</tr>
<tr>
  <td rowspan="1" colspan="1" style="text-align:center;">gmxpbsa</td>
  <td rowspan="1" colspan="1" style="text-align:center;">bash, perl</td>
  <td rowspan="1" colspan="1" style="text-align:center;">无须</td>
  <td rowspan="1" colspan="1" style="text-align:center;">Win/Lin</td>
  <td rowspan="1" colspan="1" style="text-align:center;">低</td>
  <td rowspan="1" colspan="1" style="text-align:center;">无</td>
  <td rowspan="1" colspan="1" style="text-align:center;">支持</td>
</tr>
<tr>
  <td rowspan="1" colspan="1" style="text-align:center;">g_mmpbsa</td>
  <td rowspan="1" colspan="1" style="text-align:center;">C/C++</td>
  <td rowspan="1" colspan="1" style="text-align:center;">需要, 也提供二进制版本</td>
  <td rowspan="1" colspan="1" style="text-align:center;">Lin</td>
  <td rowspan="1" colspan="1" style="text-align:center;">高</td>
  <td rowspan="1" colspan="1" style="text-align:center;">有</td>
  <td rowspan="1" colspan="1" style="text-align:center;">不支持</td>
</tr>
</table>
