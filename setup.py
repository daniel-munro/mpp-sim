from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='mppsim',
    author='Daniel Munro',
    author_email='dan@dmun.ro',
    version='1.0',
    description='Multi-parent population simulator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/daniel-munro/mpp-sim',
    packages=['mppsim'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'progress',
        'matplotlib'
    ],
    extras_require={
        'tests': [
        ],
        'docs': [
            'sphinx',
            'sphinx_rtd_theme',
            'sphinx_autodoc_typehints',
        ]
    }
)
