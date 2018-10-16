# LDPC码误码率的蒙特卡罗仿真
代码内包含基本的LDPC码编码器，译码器实现方式，以及BPSK通信系统的基本仿真原理（AWGN信道）。这被用来实现LDPC码误码率的仿真，同时也是我托管在GitHub上的第一份代码，用于学习GitHub的使用。
代码的更多说明可参见我的[博客](http://www.cnblogs.com/sea-wind2)。
## 1. 运行环境和使用方法
  运行环境：MATLAB 2014a
  使用方法：设置相应参数，运行LDPC_Simulation，可自动生成并绘制校验矩阵、生成矩阵。同时Monte Carlo仿真将运行，并将结果保存在文本文档内。Monte Carlo 仿真会消耗大量时间。（效率较低不建议用做误码率仿真）
## 2. 代码说明
  + LDCP_Simulation：脚本
  + ccsdscheckmatrix：校验矩阵构造，1、2分别对应不同文档
  + ccsdsgeneratematrix：生成矩阵构造，1、2分别对应不同文档
  + ldpcdecoderXX：对应三种不同的译码算法

## 3. 结果示例
![误码率曲线](./BER.png)
