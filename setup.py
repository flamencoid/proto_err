try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Error Correct FrameWork',
    'author': 'Phelim Bradley',
    'author_email': 'phelim.bradley@well.ox.ac.uk',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['NAME'],
    'scripts': [],
    'name': 'onec'
}

setup(**config)
