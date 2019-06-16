# 数据结构与算法大作业



**在开始写代码之前，请确认如下内容**



### 需要实现的功能

- 有一个良好的图形用户界面
    - 实现伪3D图像的显示（3D $\to$ 2D）
    - 3D图像上实现旋转，平移，缩放等操作，并可以有快捷键。（比如鼠标拖拽，滚轮等）
    - 可以计算时及时地显示Log（后台计算程序调用自定义的tk事件实现）
    - 输入由类似于Chemdraw的方法实现（个人不觉得比三维的输入简单）
    - 关于键长、键角、二面角的设定，可以由二级窗口实现，同时考虑如何实现选择。
- 实现MM（以及其他的可能（不存在的）的）算法
    - 先由python实现，尽量把涉及大量计算的部分用方法独立出来。方便可能有的C或者Fortran实现
- 实现较好的文件输入输出
    - 读入Gaussian的输入文件
    - 自行实现的输入文件及储存
    - 实现计算Log （输出文件）的储存，并且生成报告
    - 实现可以输出当前的3D图形作为图像文件



### 文件格式

参见Gaussian09的输入文件格式



### 类名及接口

暂无



```python
	if self._times == 0:  
	    self._potential = newpotential  
	else:  
	    if newpotential < self._potential:  
	        self._step *= 1.2  
	        self._potential = newpotential  
	        self.sites = self._newsites[:]  
	    else:  
	        self._step *= 0.6  
```

