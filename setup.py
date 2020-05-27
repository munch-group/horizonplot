
import setuptools

setuptools.setup(name='horizonplot',
      version='1.0.6',
      author='Kasper Munch',
      description='Generates horizon plots.',
      # long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/kaspermunch/horizonplot',
      packages=setuptools.find_packages(),
      python_requires='>=3.6',
      install_requires=[
      'pandas>=1.0',
      'numpy>=1.1',
      'seaborn>=0.1',
      'matplotlib>=3.0',
      ])
