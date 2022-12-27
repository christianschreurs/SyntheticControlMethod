from setuptools import setup
import os


#Get version
here = os.path.abspath(os.path.dirname(__file__))
_version = {}
_version_path = os.path.join(here, 'SyntheticControlMethod', '__version__.py')
with open(_version_path, 'r', 'utf-8') as f:
    exec(f.read(), _version)
    
setup(
    name='SyntheticControlMethod',
    version=_version['__version__'],
    author='Christian Schreurs',
    author_email='christian-schreurs@live.nl',
    url='https://github.com/OscarEngelbrektson/SyntheticControlMethods',
    download_url='https://github.com/christianschreurs/synthetic-control-method',
    description= "A Python package for Synthetic Control Methodology",
)
