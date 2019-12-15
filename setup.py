from setuptools import setup

setup(name='get_pbmr',
    version='0.2.4',
    description='Getting per base mutation rate from SAM and BAM files.',
    author='Semar Petrus',
    packages=['get_pbmr'],
    entry_points = {
        'console_scripts': [
            'get_pbmr=get_pbmr.get_pbmr:gen_mutation_rate_graph']},
    install_requires=['matplotlib', 'pandas', 'argparse', 'pysam'])
