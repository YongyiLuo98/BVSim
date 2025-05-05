from setuptools import setup, find_packages
from pathlib import Path

# 读取README内容作为长描述
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="BVSim",
    version="1.0.0",  # 请更新为您的版本号
    author="Yongyi LUO",
    author_email="yongyiluo98@gmail.com",
    description="Biological Variant Simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    # 包配置
    packages=find_packages(include=["main", "main.*"]),
    package_dir={"": "."},
    include_package_data=True,
    
    # 依赖项 (从您的conda环境转换而来)
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.84",
        "pysam>=0.22.1",
        "numpy>=2.1.1",
        "pandas>=2.2.2",
        "scipy>=1.14.1",
        "matplotlib>=3.9.2",
        "seaborn>=0.13.2",
        "statsmodels>=0.14.3",
        "psutil>=6.0.0",
    ],
    
    # 入口点
    entry_points={
        "console_scripts": [
            "bvsim=main.main:main",  # 创建命令行命令
        ],
    },
    
    # 分类信息
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
    ],
    
    # 项目URL
    project_urls={
        "Source": "https://github.com/YongyiLuo98/BVSim",
        "Bug Reports": "https://github.com/YongyiLuo98/BVSim/issues",
    },
)