from setuptools import setup

setup(
    name='raw_assemble_gene',
    version='1.0.0',
    description='Assemble a gene of interest from sequencing reads without assembling the full genome.',
    author='Bioinformatica-Inbiotec',
    author_email='tu@email.com',
    url='https://github.com/Bioinformatica-Inbiotec/raw_assemble_gene',
    py_modules=['assemble_gene'],
    install_requires=[
        'biopython',
        'blast',
    ],
    entry_points={
        'console_scripts': [
            'raw_assemble_gene = assemble_gene:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
)
