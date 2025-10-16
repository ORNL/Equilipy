from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os
import glob
import sys
import subprocess
from numpy import get_include

class f2py_Extension(Extension):
    def __init__(self, PackageName, srcdir):
        Extension.__init__(self, PackageName, sources=[])
        paths = sorted([d for d in glob.glob(f'{srcdir[0]}/*/') if os.path.isdir(d)])
        
        sourcedirs = []
        for path in paths:
            filenames = glob.glob(f'{path}/*.f*')
            sourcedirs.extend(filenames)
            
        self.sourcedirs = [os.path.abspath(sourcedir) for sourcedir in sourcedirs]
        self.dirs = list(sourcedirs)
        self.include_dirs = [get_include()]
        self.packages = find_packages(where=PackageName)

class f2py_Build(build_ext):
    slist = ''

    def run(self):
        for ext in self.extensions:
            self.compile_fortran(ext)

        f2py_command = [
            sys.executable, '-m', 'numpy.f2py', '-c',
            '-m', 'equilifort', '--backend', 'meson'
        ]
        f2py_command.extend(self.slist.split())

        print(f"Running command: {' '.join(f2py_command)}")
        
        result = subprocess.run(f2py_command, cwd='equilipy', capture_output=True, text=True)

        print("f2py stdout:\n", result.stdout)
        print("f2py stderr:\n", result.stderr)

        if result.returncode != 0:
            raise Exception("f2py compilation failed.")

        lib_name_base = self.get_ext_filename('equilifort').split('/')[-1]
        
        built_lib_path = None
        for root, dirs, files in os.walk('equilipy'):
            if lib_name_base in files:
                built_lib_path = os.path.join(root, lib_name_base)
                print(f"Found compiled library at: {built_lib_path}")
                break
        
        if not built_lib_path:
            raise FileNotFoundError(f"Could not find {lib_name_base} after f2py build.")

        build_lib = os.path.abspath(self.build_lib)
        final_lib_path = os.path.join(build_lib, 'equilipy', lib_name_base)
        
        os.makedirs(os.path.dirname(final_lib_path), exist_ok=True)
        print(f"Copying {built_lib_path} to {final_lib_path}")
        self.copy_file(built_lib_path, final_lib_path)
        
        self.copy_file(built_lib_path, os.path.join('equilipy', lib_name_base))
        
        for f in glob.glob('equilipy/*.mod'):
            os.remove(f)
        
    def compile_fortran(self, ext):
        for file in ext.dirs:
            rel_path = os.path.relpath(file, 'equilipy')
            self.slist += f'{rel_path} '
            
setup(
    packages=find_packages(),
    ext_modules=[f2py_Extension('equilipy', ['equilipy/fsrc'])],
    cmdclass=dict(build_ext=f2py_Build),
)
