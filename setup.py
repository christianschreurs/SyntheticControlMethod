from setuptools import setup
import os

install_requires = [
        'numpy>=1.17',
        'scipy>=1.4.1',
        'pandas>=1.1.2',
        'matplotlib>=2.2.3',
    ]

"""
here = os.path.abspath(os.path.dirname(__file__))
_version = {}
_version_path = os.path.join(here, 'SyntheticControlMethod', '__version__.py')
with open(_version_path, 'r', 'utf-8') as f:
    exec(f.read(), _version)
"""
    
setup(
    name='SyntheticControlMethod',
    author='Christian Schreurs',
    author_email='christian-schreurs@live.nl',
    url='https://github.com/OscarEngelbrektson/SyntheticControlMethods',
    download_url='https://github.com/christianschreurs/synthetic-control-method',
    description= "A Python package for Synthetic Control Methodology",
    install_requires=install_requires
)
