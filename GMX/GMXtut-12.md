---
 layout: post
 title: GROMACS教程：Xmgrace学习笔记
 categories:
 - 科
 tags:
 - GMX
---

* toc
{:toc}


<ul class="incremental">
<li>整理: 李卫星</li>
</ul>

<p>对于用惯Linux的人, 在Linux下画图有点不方便, gnuplot虽然挺不错, 但是不能直接在图上进行操作, 这是它的不足吧. 为此, 简单介绍一下xmgrace.</p>

## 安装

<p>xmgrace可以直接用命令来安装：</p>

<pre><code>sudo apt-get install xmgrace
</code></pre>

<p>也可以下载软件包自己编译. 官网<a href="http://plasma-gate.weizmann.ac.il/Grace/">http://plasma-gate.weizmann.ac.il/Grace/</a>.</p>

<p><a href="xmgrace.zip">这里</a>有一个整理好的压缩包, 里面包含了xmgrace&#8211;5.1.25的源代码以及一些资料.</p>

<p>grace和xmgrace差不多, 具体区别自己看网上.</p>

## 一些资源

<p>几个网址：</p>

<ol class="incremental">
<li><a href="http://plasma-gate.weizmann.ac.il/Grace/phpbb/index.php?sid=df4c3ea5b571195761e69142e08939d3">Discussions讨论区</a></li>
<li><a href="http://mintaka.sdsu.edu/reu/grace.tutorial.html#toc1">An Xmgrace Tutorial</a></li>
<li><a href="http://plasma-gate.weizmann.ac.il/Grace/doc/Tutorial.html">Grace Tutorials</a>, 个人觉得这个是不错的</li>
<li><a href="http://math.nyu.edu/aml/software/xmgrace.html">Grace/Xmgrace (Xmgr)</a></li>
<li><a href="http://plasma-gate.weizmann.ac.il/Grace/doc/UsersGuide.html">Grace User&#8217;s Guide</a></li>
</ol>

<p>下载上面的grace&#8211;5.1.25后, Grace Tutorial教程实例文件在<code>grace-5.1.25/doc</code>下面（Tutorial.pdf也在）. 还有一些做好的图例在<code>grace-5.1.25/examples</code>, 从那里可以看到xmgrace也是可以的. 下面是随便挑的几张, 图片的的字体都是用xmgrace加上去的.</p>

<figure>
<img src="/GMX/xmgrace_1.png" alt="" />
</figure>

<figure>
<img src="/GMX/xmgrace_2.png" alt="" />
</figure>

<p>本文档主要是学习Grace Tutorials的总结, 如有不对, 请指正.</p>

## 基本操作

### 1. 简单画图：使用数据文件

<p>可以在工作目录下打开终端, 输入<code>xmgrace file</code>.</p>

<p>如果数据文件里面是多列数据, 即X Y1 Y2&#8230;, 可以输入<code>xmgrace -nxy  file</code>. 加上<code>-nxy</code>就可以是多条曲线了.</p>

<p>如果要对坐标轴取对数, 可以添加选项 <code>-log x|y|xy</code>, 如：<code>xmgrace -log x -nxy file</code>. （好像要放前面, 放后面会出错！）</p>

### 2. 简单绘图：导入文件

<p>也可以在终端只输入<code>xmgrace</code>, 不带文件名, 打开空白的xmgrace再自己导入. 点击菜单<code>Data | Import | ASCII</code>在<code>files</code>下拉框选择文件, 点一下<code>OK</code>就可以了.</p>

<p>注意, 好像只有扩展名为<code>.dat</code>的才显示, 而<code>.xvg</code>格式的不显示, 需要自己改一下扩展名. 还有就是, 如果文件路径有中文的话显示会出现乱码, 但仍然可以选择, 只要你确定哪个是你的文件.</p>

### 3. 创建数据表格

<p>也可以自己用表格创建绘图数据. 在终端输入<code>xmgrace</code>, 打开空白的xmgrace. 点击菜单<code>Main | Edit/Data_sets...</code>窗口, 在上面的空白框中, 一直右击按住出现菜单移动（右击还不要放开）, 选<code>creat new</code>选<code>in spreadsheet</code>（这里还有另外三个, 如编写公式输入, 你可以自己摸索或看手册）打开表格, 就可以输入x y点数值了.</p>

### 4. 保存图像

<p>画好了图片就保存, 点击<code>File | print_setup</code>, 在<code>device</code>选格式, 我常选<code>png</code>, 对不起的是好像没有<code>tiff</code>格式可选. 输入想要保存的名字, 点<code>aceept</code>退出窗口, 按一下快捷键<code>ctrl+p</code>就保存了. 其实也可以用<code>File | print</code>, 但最好记住快捷键. 看看目录内是否有图片文件.</p>

<p>xmgrace不像origin那样有什么工程文件, 它的设置不会自动保存, 所以你上一次对数据的操作是不会保存的. 这很麻烦. 所以如果我们以后还要对文件进行操作（如坐标轴, 曲线粗细不够的等问题）, 很麻烦, 其他的设置要重新做一遍. 为避免这些麻烦, 就得保存设置. 点击菜单<code>File | save</code>或<code>File | save as</code>. <code>File | save</code>直接把设置加入到你的数据文件开头, 自己操作完可以去看看. <code>File | save as</code>就是另存一份, 原数据文件开头添加设置. 如果开头保存有设置的话, 直接右击用xmgrace查看, 就可以看见你原来设置的图片. 你修改后, 重新<code>File | save</code>. 这样也有一些好处, 如果多个文件都是一样设置, 做好一个文件, 把设置保存到数据文件开头, 然后把它们复制到其他文件, 就搞定了, 然后再用xmgrace输出图片.</p>

<p>你可以去<code>grace-5.1.25/exemple</code>打开各个文件学习一下各种图是如何进行设置的.</p>

### 5. 选择数据列

<p>第一点说过如何使用<code>-nxy</code>选项绘制多列数据, 但如果只想用第1和第3列数据画图, 咋办? 可以使用<code>xmgrace -block file.dat -bxy 1:3</code>. 其中<code>-bxy</code>就是选1和3两列.
也可以在<code>-nxy</code>打开多条曲线后, 点击图上的曲线, 出现<code>set appearance</code>窗口, 在<code>select set</code>下, 选中不要的数据列, 长按右击移动选<code>hide</code>就可以了（数据列多的话, 就麻烦一些）;</p>

<p>说到这里, 说一下对文件的数据处理, 看示例：</p>

<pre><code>xmgrace -nosafe -nxy box.xvg -pexec &quot;s0.y=(s0.y*s1.y)/64&quot; -pexec &quot;kill s1&quot; -pexec &quot;autoscale&quot;
</code></pre>

<p>我的文件<code>box.xvg</code>中有三列数据, 时间帧、盒子x大小、盒子y大小, 计算膜表面的APL(area prea lipid)是用x*y/64. 用xmgrace怎么计算呢? 像上面那样用<code>-pexec</code>选项输入参数命令就可以了. 输出的图像就是x轴是时间, y轴是经过处理的得到的APL. 更具体的请看手册教程.</p>

<p>还是再说一些多列数据的情况. <code>data | Import/ASCII</code>选好dat文件后在中间那里看到<code>load as</code>了吧, 点一下选<code>block data</code>, 后点<code>ok</code>跳出来框框, <code>x from column</code>自己选<code>y from column</code>, 选好就可以. 如果都选1, 会是什么图像, 猜一下, xy都是一样的值, 当然是45度斜率的直线y=x. 那里可以选偏差条, 就是xydy, 自己摸索了.</p>

### 6. 图像设置

<p>数据选择操作说的差不多, 说一下图片的外观设置（也就是具体的坐标轴文字大小、间距、范围, 图片标题, 线条粗、细颜色等）. 我们可以双击相应的地方（和origin类似, 自己体会）也可以菜单里选择<code>plot | plot appearance</code>.</p>

<p>说一下一些基本的操作：</p>

<ul class="incremental">
<li><code>Legends</code>: 图例, 就是多条曲线时小框用不同颜色标记出来</li>
<li><code>frame</code>: 画出的图的框框</li>
<li><code>tick label</code>: 就是坐标的下面的坐标间隔标度1 2 3 4 5 6</li>
<li><code>tick mark</code>: 就是坐标轴的突起的刻度, 如 1 头上对应轴有个突起</li>
<li><code>leg.box</code>: 在标题窗口里面. 可以调整图例的位置, 里面的<code>location</code>设置就可以了</li>
<li><code>axis placement</code>: 用处在 y 轴 坐标在&#8211;1.5<sub>1</sub>.5 范围时, 坐标轴若在y=&#8211;1.5 处 就去<code>zreo axis</code>点一下就可以了.</li>
</ul>

<p>图像的四个角可以拉大拉小的, 如果不合适自己调.</p>

<p>如果有多条曲线, 出现了图例与曲线重叠, 可以试试这样操作：拖动图例的方法, <code>Ctrl+L</code>单击, 就可以使箭头变手形, 拖动图例了. 如还不行, 双击再试试看可以了吗. 也可以在<code>appearance</code>更改.</p>

<p>记住<code>set appearance</code>要在窗口里面的<code>select set</code>位置选相应的<code>s0 s1....</code>, <code>main</code>标签下的<code>Legend</code>就是图例, 你可以给每条曲线标记不同名字以便区分.</p>

<p>还有需要说一下的画布左边的按钮, 就说一个吧, <code>AS</code>表示恢复（图片设置的外观不变的, 好）. 所以你按钮按错了, 点它就行了.</p>

<p>如果用不同坐标刻度值时, 也就是跳跃很大, 如μs, ms, s, 在<code>Axes</code>（双击坐标轴就出来了）的窗口的<code>tick properties</code>框里面选<code>Format</code>, 里面好几种格式, 选<code>Compute(K, M, G....)</code>或<code>Engnieering</code>, 这样就可以一个坐标轴多个单位了. 不用显示那么长的数字串.</p>

### 7. 多图并列

<p>有时候要一个画布放多张图, 打开菜单<code>Edit | Arrange_graphs</code>对话框, 填写需要几rows几columns, 点应用, 就出来四个小框了（假设2*2格式）. 具体外观, 可以在刚才的框框里继续调节, 点击各个小框点击, 像前面的第2点那样选好文件进来就可以. 具体其他调节参看前面的6步骤操作.</p>

<p>也可以使用命令输入, 加上指令就可以了：</p>

<pre><code>-pexec &quot;arrange (2, 2, .1, .1, .1, ON, ON, ON)&quot;
</code></pre>

<p>前面的两个2代表2*2排列, 加上<code>-graph</code>更好, 代表的是后面接的文件放在哪个位置, 从0开始数：</p>

<pre><code>xmgrace -graph 2 10_rdf.xvg -graph 1 11.dat -graph 0 rdf01.dat -graph 3 rdf.dat -pexec &quot;arrange(2,2,.1,.1,.1,ON,ON,ON )&quot;
</code></pre>

### 8. 内嵌图形

<p>说实话, xmgrace做内嵌图像, 我也没懂, 希望会的人能介绍一下. 参考方法：</p>

<p>打开菜单<code>Edit | overlay_graphs</code>, 点中在<code>overlay graph</code>的<code>G0</code>, 然后右击, 选<code>creat new</code>, 出现了另一个空数据集. 然后<code>overlay graph</code>和<code>on to</code>, 反正个数据集就可以, 后点同意, 就两个图的, 当然看起来是一个, 自己点一下四周点移动, 就可以看见两个的, 只是重叠而已. 然后一个一个添加数据了. 刚才说了内嵌图形还记得吗, 这里也可以把其中一个图拉小移动作为一个内嵌图, 再修改一下需要显示的坐标范围. 具体靠自己去实践了.</p>

### 9. 双y轴曲线

<p>有时候我们需要同一图像两条不同y轴的曲线, 如例子：坐标x表示时间, 左边y轴代表压力, 右边y轴代表是温度. 注意调整坐标刻度, 在<code>Axes</code>窗口<code>tick label</code>和<code>tick mark</code>的<code>draw  on</code>, 一个选正常<code>normal side</code>, 一个选<code>oppsite side</code>相反, 就可以了. 对边的刻度就不要显示了.</p>

### 10. 数据拟合

<p>选项在<code>data | tramsformation</code>.</p>

<p><strong>线性拟合</strong></p>

<p>打开图像后</p>

<ol class="incremental">
<li>点击菜单<code>DATA | TRANSFORMATION | REGRESSION</code></li>
<li>选择<code>SET</code></li>
<li>选择<code>LINEAR FIT</code></li>
<li>按下<code>ACCEPT</code>n</li>
<li>出来一个窗口显示拟合的直线方程表达式了. 如果需要斜率, 就自己记下方程. Save the information about the slope and intercept that appear in a blue console as slope.dat.</li>
</ol>

<p><strong>二次方程</strong></p>

<p>点击<code>DATA | TRANSFORMATION | INTERPOLATION/SPLINE</code></p>

<p>选择<code>SET | METHOD | CUBIC SPLINE</code>, <code>START</code>设为1, <code>STOP</code>设为10, <code>LENGTH</code>设为500或1000, 然后点击`ACCEPT.</p>

<p><strong>非线性拟合</strong></p>

<p>输入方程形式, 还要输入参数初始值.</p>

## 附录: 图像文本设置

<p>如果像gnuplot文件那样修改设置的话, 我们直接修改文本就可以了. 可以参考<code>grace-5.1.25/exemple</code>下的示例文件去学习, 吃透了就是xmgrace的高手了. 下面是其中部分参数的意义.</p>

<ul class="incremental">
<li>@ title &quot;&quot; 标题</li>
<li>@ title font 0 标题字体格式</li>
<li>@ title size 1.500000 标题字体大小</li>
<li>@ xaxis bar linewidth 3.0 坐标轴粗细</li>
<li>@ xaxis label &#8220;E （KJ/mol）&#8221; x坐标轴表示的物理量</li>
<li>@ xaxis tick major 100 单位刻度100</li>
<li>@ yaxis label &quot;&quot; y坐标轴表示的物理量</li>
<li>@ yaxis tick major 2 每隔2画一个标度, 也就是精度为1</li>
<li>@ s5 line linewidth 3.0 第6条曲线的的粗细</li>
<li>@ xaxis tick major 轴上的标度</li>
<li>@ xaxis ticklabel char size 2.000000 轴上的标度 显示大小</li>
<li>@ s3 hidden false 是否隐藏曲线</li>
<li>@ s3 legend &#8220;Interface&#8221; 第四条曲线的图例的名称标示, 就是表示哪条曲线</li>
<li>@ s3 symbol 5 曲线的点的符号, 如圆圈, 星号, 倒三角</li>
<li>@ s3 line type 5 曲线变成长虚线, 1是直线</li>
<li>@ page size 792, 612 白色画布的大小</li>
</ul>

<p>颜色代码的意义</p>

<ul class="incremental">
<li>@map color 0 to (255, 255, 255), &#8220;white&#8221;</li>
<li>@map color 1 to (0, 0, 0), &#8220;black&#8221;</li>
<li>@map color 2 to (255, 0, 0), &#8220;red&#8221;</li>
<li>@map color 3 to (0, 255, 0), &#8220;green&#8221;</li>
<li>@map color 4 to (0, 0, 255), &#8220;blue&#8221;</li>
<li>@map color 5 to (255, 255, 0), &#8220;yellow&#8221;</li>
<li>@map color 6 to (188, 143, 143), &#8220;brown&#8221;</li>
<li>@map color 7 to (220, 220, 220), &#8220;grey&#8221;</li>
<li>@map color 8 to (148, 0, 211), &#8220;violet&#8221;</li>
<li>@map color 9 to (0, 255, 255), &#8220;cyan&#8221;</li>
<li>@map color 10 to (255, 0, 255), &#8220;magenta&#8221;</li>
<li>@map color 11 to (255, 165, 0), &#8220;orange&#8221;</li>
<li>@map color 12 to (114, 33, 188), &#8220;indigo&#8221;</li>
<li>@map color 13 to (103, 7, 72), &#8220;maroon&#8221;</li>
<li>@map color 14 to (64, 224, 208), &#8220;turquoise&#8221;</li>
<li>@map color 15 to (0, 139, 0), &#8220;green4&#8221;</li>
</ul>

<p>上面的这些在<code>.xvg</code>文件开头写出来给我们参考. 让我们明白这些数字代表什么颜色. 还有下面的字体格式部分</p>

<p>字体代码的意义：</p>

<ul class="incremental">
<li>@map font 8 to &#8220;Courier&#8221;, &#8220;Courier&#8221;</li>
<li>@map font 10 to &#8220;Courier-Bold&#8221;, &#8220;Courier-Bold&#8221;</li>
<li>@map font 11 to &#8220;Courier-BoldOblique&#8221;, &#8220;Courier-BoldOblique&#8221;</li>
<li>@map font 9 to &#8220;Courier-Oblique&#8221;, &#8220;Courier-Oblique&#8221;</li>
<li>@map font 14 to &#8220;Courier-Regular&#8221;, &#8220;Courier-Regular&#8221;</li>
<li>@map font 15 to &#8220;Dingbats-Regular&#8221;, &#8220;Dingbats-Regular&#8221;</li>
<li>@map font 4 to &#8220;Helvetica&#8221;, &#8220;Helvetica&#8221;</li>
<li>@map font 6 to &#8220;Helvetica-Bold&#8221;, &#8220;Helvetica-Bold&#8221;</li>
<li>@map font 7 to &#8220;Helvetica-BoldOblique&#8221;, &#8220;Helvetica-BoldOblique&#8221;</li>
<li>@map font 5 to &#8220;Helvetica-Oblique&#8221;, &#8220;Helvetica-Oblique&#8221;</li>
<li>@map font 20 to &#8220;NimbusMonoL-Bold&#8221;, &#8220;NimbusMonoL-Bold&#8221;</li>
<li>@map font 21 to &#8220;NimbusMonoL-BoldOblique&#8221;, &#8220;NimbusMonoL-BoldOblique&#8221;</li>
<li>@map font 22 to &#8220;NimbusMonoL-Regular&#8221;, &#8220;NimbusMonoL-Regular&#8221;</li>
<li>@map font 23 to &#8220;NimbusMonoL-RegularOblique&#8221;, &#8220;NimbusMonoL-RegularOblique&#8221;</li>
<li>@map font 24 to &#8220;NimbusRomanNo9L-Medium&#8221;, &#8220;NimbusRomanNo9L-Medium&#8221;</li>
<li>@map font 25 to &#8220;NimbusRomanNo9L-MediumItalic&#8221;, &#8220;NimbusRomanNo9L-MediumItalic&#8221;</li>
<li>@map font 26 to &#8220;NimbusRomanNo9L-Regular&#8221;, &#8220;NimbusRomanNo9L-Regular&#8221;</li>
<li>@map font 27 to &#8220;NimbusRomanNo9L-RegularItalic&#8221;, &#8220;NimbusRomanNo9L-RegularItalic&#8221;</li>
<li>@map font 28 to &#8220;NimbusSansL-Bold&#8221;, &#8220;NimbusSansL-Bold&#8221;</li>
<li>@map font 29 to &#8220;NimbusSansL-BoldCondensed&#8221;, &#8220;NimbusSansL-BoldCondensed&#8221;</li>
<li>@map font 30 to &#8220;NimbusSansL-BoldCondensedItalic&#8221;, &#8220;NimbusSansL-BoldCondensedItalic&#8221;</li>
<li>@map font 31 to &#8220;NimbusSansL-BoldItalic&#8221;, &#8220;NimbusSansL-BoldItalic&#8221;</li>
<li>@map font 32 to &#8220;NimbusSansL-Regular&#8221;, &#8220;NimbusSansL-Regular&#8221;</li>
<li>@map font 33 to &#8220;NimbusSansL-RegularCondensed&#8221;, &#8220;NimbusSansL-RegularCondensed&#8221;</li>
<li>@map font 34 to &#8220;NimbusSansL-RegularCondensedItalic&#8221;, &#8220;NimbusSansL-RegularCondensedItalic&#8221;</li>
<li>@map font 35 to &#8220;NimbusSansL-RegularItalic&#8221;, &#8220;NimbusSansL-RegularItalic&#8221;</li>
<li>@map font 36 to &#8220;StandardSymbolsL-Regular&#8221;, &#8220;StandardSymbolsL-Regular&#8221;</li>
<li>@map font 12 to &#8220;Symbol&#8221;, &#8220;Symbol&#8221;</li>
<li>@map font 38 to &#8220;Symbol-Regular&#8221;, &#8220;Symbol-Regular&#8221;</li>
<li>@map font 2 to &#8220;Times-Bold&#8221;, &#8220;Times-Bold&#8221;</li>
<li>@map font 3 to &#8220;Times-BoldItalic&#8221;, &#8220;Times-BoldItalic&#8221;</li>
<li>@map font 1 to &#8220;Times-Italic&#8221;, &#8220;Times-Italic&#8221;</li>
<li>@map font 0 to &#8220;Times-Roman&#8221;, &#8220;Times-Roman&#8221;</li>
<li>@map font 43 to &#8220;URWBookmanL-DemiBold&#8221;, &#8220;URWBookmanL-DemiBold&#8221;</li>
<li>@map font 44 to &#8220;URWBookmanL-DemiBoldItalic&#8221;, &#8220;URWBookmanL-DemiBoldItalic&#8221;</li>
<li>@map font 45 to &#8220;URWBookmanL-Light&#8221;, &#8220;URWBookmanL-Light&#8221;</li>
<li>@map font 46 to &#8220;URWBookmanL-LightItalic&#8221;, &#8220;URWBookmanL-LightItalic&#8221;</li>
<li>@map font 47 to &#8220;URWChanceryL-MediumItalic&#8221;, &#8220;URWChanceryL-MediumItalic&#8221;</li>
<li>@map font 48 to &#8220;URWGothicL-Book&#8221;, &#8220;URWGothicL-Book&#8221;</li>
<li>@map font 49 to &#8220;URWGothicL-BookOblique&#8221;, &#8220;URWGothicL-BookOblique&#8221;</li>
<li>@map font 50 to &#8220;URWGothicL-Demi&#8221;, &#8220;URWGothicL-Demi&#8221;</li>
<li>@map font 51 to &#8220;URWGothicL-DemiOblique&#8221;, &#8220;URWGothicL-DemiOblique&#8221;</li>
<li>@map font 52 to &#8220;URWPalladioL-Bold&#8221;, &#8220;URWPalladioL-Bold&#8221;</li>
<li>@map font 53 to &#8220;URWPalladioL-BoldItalic&#8221;, &#8220;URWPalladioL-BoldItalic&#8221;</li>
<li>@map font 54 to &#8220;URWPalladioL-Italic&#8221;, &#8220;URWPalladioL-Italic&#8221;</li>
<li>@map font 55 to &#8220;URWPalladioL-Roman&#8221;, &#8220;URWPalladioL-Roman&#8221;</li>
<li>@map font 56 to &#8220;Utopia-Bold&#8221;, &#8220;Utopia-Bold&#8221;</li>
<li>@map font 57 to &#8220;Utopia-BoldItalic&#8221;, &#8220;Utopia-BoldItalic&#8221;</li>
<li>@map font 58 to &#8220;Utopia-Italic&#8221;, &#8220;Utopia-Italic&#8221;</li>
<li>@map font 59 to &#8220;Utopia-Regular&#8221;, &#8220;Utopia-Regular&#8221;</li>
<li>@map font 13 to &#8220;ZapfDingbats&#8221;, &#8220;ZapfDingbats&#8221;</li>
</ul>

<p>最后还是说一句, 有些方法不懂话, 看看<code>grace-5.1.25/examples</code>下面的例子, 右击用grace打开, 看看人家的图的设置, 有你想要的效果吗, 仔细琢磨一下人家怎么设置的. 这是最好的学习方法, 比看教程好很多, 当然这是个人认为.</p>