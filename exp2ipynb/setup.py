import setuptools

setuptools.setup(
    name='Exp2Ipynb',
    version='0.0.5',
    Author='Ulf Liebal',
    author_email='ulf.liebal@rwth-aachen.de',
    description='Exp2Ipynb supports the analysis of gene expression libraries.',
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy>=1.18.0',
        'pandas',
        'seaborn',
        'matplotlib>=3.1.0',
        'joblib>=0.14.0',
        'scipy>=1.4.0',
        'pylab-sdk>=1.1.0',
        'openpyxl>=3.0.0',
        'scikit-learn>=0.22.0',
        'logomaker>=0.8',
        'deap>=1.3.0',
        'pydot>=1.4.1'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GPL3',
        'Operating System :: OS Independent'
    ],
    python_requires='>=3.7'
)
