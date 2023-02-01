from distutils.core import setup, Extension
import numpy

def main():

    setup(name="ssa",
          version="1.0.0",
          description="C library functions for SSA",
          author="Clayton Seitz",
          author_email="cwseitz@uchicago.edu",
          ext_modules=[Extension("ssa._ssa", ["ssa/_ssa/ssa.c"],
                       include_dirs = [numpy.get_include(), '/usr/include/gsl'],
                       library_dirs = ['/usr/lib/x86_64-linux-gnu'],
                       libraries=['m', 'gsl', 'gslcblas'])])


if __name__ == "__main__":
    main()
