from setuptools import setup

config = {
    'description': 'Clonify',
    'author': 'Bryan Briney',
    'url': 'www.github.com/briney/clonify/',
    # 'download_url': 'www.github.com/briney/abcloud/',
    'author_email': 'briney@scripps.edu',
    'version': '0.1.0',
    'install_requires': ['nose',
                         # 'abtools',
                         ],
    'packages': ['clonify'],
    'scripts': ['bin/clonify'],
    'name': 'clonify',
    'include_package_data': True
}

setup(**config)
