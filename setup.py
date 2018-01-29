from distutils.core import setup

setup(
    name='cmm',
    version='0.99',
    author="Haohan Wang",
    author_email='haohanw@cs.cmu.edu',
    url="https://github.com/HaohanWang/CMM",
    description="Joint Genetic Analysis of Complex Disorders from Independently Collected Data Sets: Application to Alzheimer's Disease and Substance Use Disorder",
    packages=['models', 'utility'],
    scripts=['cmm.py'],
)
