<!-- Author : Changxing Su-->

<div align="center">

# MicroCAT

</div>

## 简介

MicroCAT是一个专为微生物组学数据分析设计的综合性计算工具箱，可以从细胞分辨率层面提取微生物信息。该工具箱整合了多种分析方法，帮助研究人员更高效地处理和解析微生物数据。

## 目录

- [简介](#简介)
- [软件安装](#软件安装)
- [使用文档](#使用文档)
- [引用](#引用)
- [联系方式](#联系方式)
- [社区贡献](#社区贡献)


## 📦 软件安装

MicroCAT可以在Python 3.8及以上版本中运行，我们提供多种安装方式：

### 1. Conda安装（推荐）
```
conda create -n MICROCAT -c bioconda microcat
```
如果您的运行环境中没有conda，我们建议您先安装conda，具体安装方法可以参考[conda的官方文档](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)。

### 2. Pip安装
使用[pip](https://pip.pypa.io/en/stable/installation/)从[PyPI](https://pypi.org/)快速安装microcat：
```
pip install microcat
```
然后安装microcat运行所需要使用的软件，或者在运行时使用'--use-conda'参数来自动构建运行环境（详见[microcat的官方文档](https://github.com/zhaofangyuan98/MicroCAT/wiki/MicroCAT-Tutorial)）。

### 3. Docker 镜像

> Docker镜像还在构建中，请耐心等待。

## 📄 说明文档

具体使用方法可以参考[microcat的官方文档](https://github.com/zhaofangyuan98/MicroCAT/wiki/MicroCAT-Tutorial)。

## 📄 引用

如果您在研究中使用了MicroCAT，请引用我们的论文（即将发布）。

```bibtex
@article{microcat,
  title={MicroCAT: Microbial Information from Cell Resolution Omics Computational Analysis Toolbox},
  author={},
  journal={},
  year={2025}
}
```

## 📧 联系方式

如有任何问题，请随时联系我们。

Email: [changxingsu42@gmail.com](mailto:changxingsu42@gmail.com) 或者 [12231359@mail.sustech.edu.cn](mailto:12231359@mail.sustech.edu.cn) 😃


## 💪社区贡献

我们欢迎您参与到microcat的开发中来，您可以通过以下方式贡献：

- 报告问题和提出建议
- 提交代码改进
- 完善文档
- 分享您的使用经验

具体贡献方法可以参考[microcat的官方文档](https://github.com/zhaofangyuan98/MicroCAT/wiki/MicroCAT-Tutorial)。

感谢所有贡献者的支持与帮助！






