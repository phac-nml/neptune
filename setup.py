from setuptools import find_packages, setup
#import neptune.Neptune

dependencies = ['drmaa', 'numpy', 'scipy', 'biopython']

setup(
    name='neptune',
    version='1.2.1',
    url='https://github.com/phac-nml/neptune.git',
    license='Apache-2.0',
    author='Eric Marinier',
    author_email='eric.marinier@phac-aspc.gc.ca',
    description='Neptune signature discovery.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'neptune = neptune.Neptune:main',
        ],
    },
)
