from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup(
  name = "pypaillier",
  ext_modules=[
    Extension("pypaillier", [ "pypaillier.pyx" ],
		include_dirs=['../libpaillier/.'], 
		library_dirs=['../libpaillier/.'],
		libraries=['gmp', "paillier"] )
    ],
  cmdclass = {'build_ext': build_ext}
)

