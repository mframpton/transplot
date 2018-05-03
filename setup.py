from setuptools import setup, find_packages

setup(
    name='transplot',
    version='1.0',
    packages=find_packages(),
    license='GNU General Public License v3.0',
    long_description=open('README.md').read(),
    author="Matthew Frampton",
    author_email="mjeframpton@gmail.com",
    description="The transplot package can be used to make multi-track plots of Next Generation Sequencing (NGS) data for gene transcripts, namely for depth of coverage and for the distribution of variants and protein domains.",
    url='https://github.com/mframpton/transplot',
    install_requires=['regex==2017.11.9','pandas==0.21.0','matplotlib==1.5.3','numpy==1.11.3'],
)