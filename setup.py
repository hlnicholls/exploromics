from setuptools import setup, find_packages 
import io

setup(name='exploromics',
	version='0.1.0',
	author='Hannah Nicholls',
	description='Multi-omic data integration and exploratory data analysis for a gene list',
	url="https://github.com/hlnicholls/exploromics",
	packages=find_packages(),
	python_requires='>=3.6',
	install_requires=['setuptools==70.0.0',
    'lxml==4.8.0',
    'matplotlib==3.6.2',
    'mygene==3.2.2',
    'numpy==1.22.3',
    'pandas==1.4.2',
    'seaborn==0.11.2'
],
	include_package_data=True,
	entry_points={'console_scripts': ['exploromics=exploromics.__main__:main',]})