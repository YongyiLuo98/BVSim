# BVSim

BVSim 是一个用于基因模拟和结构变异的 Python 包。

## 安装

你可以使用 pip 来安装 BVSim：

\`\`\`bash
pip install BVSim
\`\`\`

## 使用

首先，你需要导入 BVSim 包：

\`\`\`python
import BVSim
\`\`\`

然后，你可以使用 BVSim 包中的函数和类来进行基因模拟和结构变异。

## 功能

- 基因模拟：BVSim 可以模拟基因序列，包括各种类型的基因和非编码区域。
- 结构变异：BVSim 可以模拟各种类型的结构变异，包括插入、删除、倒位和转位。

## 贡献

如果你对 BVSim 有任何建议或问题，欢迎提交 issue 或 pull request。

## 许可

BVSim 使用 MIT 许可证，你可以在 LICENSE 文件中查看详细信息。

# 安装和使用指南

## 安装

1. 克隆此仓库：`git clone https://github.com/yourusername/BVSim.git`
2. 进入仓库目录：`cd BVSim`
3. 安装必要的依赖：`pip install -r requirements.txt`

## 使用

运行BVSim：

```bash
python -m BVSim -ref '/path/to/reference.fasta' -save '/path/to/save/directory/' -seed 1 -rep 5 -write -snp 2000
```
