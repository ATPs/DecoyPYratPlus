from setuptools import setup

version = '2.0.1'

setup(
    name = 'DecoyPYratPlus',
    version = version,
    packages = ['DecoyPYratPlus'],
    description = 'Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectrometery Analyses',
    author = 'James Wright, Xiaolong Cao',
    author_email = 'atps@outlook.com',
    url = 'https://github.com/ATPs/DecoyPYratPlus',
    keywords = ['database searching',' fdr',' python',' sequence database',' shotgun proteomics',' target-decoy'],
    install_requires=['numpy',
    ],
    license = 'MIT',
    entry_points={
        'console_scripts': [
            'DecoyPYratPlus = decoypyrat.decoyPYratPlus:main'
        ]
    }
)
