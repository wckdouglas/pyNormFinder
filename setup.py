import sys
from setuptools import setup, find_packages
import glob

setup(
    name='NormFinder',
    author='Douglas Wu',
    author_email='douglas.wu@qiagen.com',
    description='Python version of NormFinder',
    license='MIT',
    packages=find_packages(),
    scripts = glob.glob('bin/*'),
    package_data={'CLCbio': ['snakefiles/*.smk']},
    install_requires=[
        'pandas',
        'numpy',
        'logging']
)
