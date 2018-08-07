from setuptools import setup

setup(name='mddspls',
      version='0.1',
      description='The multi data driven sparse pls package',
      url='http://github.com/hlorenzo/mddspls',
      author='Hadrien Lorenzo',
      author_email='hadrien.lorenzo.2015@gmail.com',
      license='MIT',
      packages=['mddspls'],
      install_requires=['numpy','pandas','warnings','sklearn',
      'multiprocessing','itertools','random'],
      zip_safe=False)