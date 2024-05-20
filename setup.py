from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os, glob, sys
from numpy import get_include

# Define the Extension for the Fortran script
class f2py_Extension(Extension):

    def __init__(self, PackageName, srcdir):
        Extension.__init__(self, PackageName, sources=[])
        # Use glob to get all directory names
        paths = sorted([d for d in glob.glob(f'{srcdir[0]}/*/') if os.path.isdir(d)])
        
        # Get paths for all fortran files
        sourcedirs = list([])
        for i,path in enumerate(paths):
            filenames = list(glob.glob(f'{path}/*.f*'))
            sourcedirs.extend(filenames)
        self.sourcedirs = [os.path.abspath(sourcedir) for sourcedir in sourcedirs]
        self.dirs = list(sourcedirs)
        self.include_dirs=[get_include()]
        self.packages=find_packages(where=PackageName)

        

class f2py_Build(build_ext):
    
    slist=''
    def run(self):
        # Compile Fortran script using f2py and create shared object file
        for ext in self.extensions:
            self.compile_fortran(ext)
            
        # Copy the shared object file to the build directory
        if sys.platform.startswith('win'):
            os.system(f'cd equilipy&& f2py -c {self.slist} -m equilifort --backend distutils --fcompiler=gnu95')
        else:
            os.system(f'cd equilipy; f2py -c {self.slist} -m equilifort --backend meson --fcompiler=gnu95')
        
        lib_name = self.get_ext_filename('equilipy/equilifort')
        
        build_lib = os.path.abspath(self.build_lib)
        lib_path = os.path.join(build_lib, lib_name)
        self.copy_file(lib_name, lib_path)
        
        # Add dll in windows platform when meson is not in use        
        if sys.platform.startswith('win'):
            dll=glob.glob('equilipy\\equilifort\\.libs\\*.dll')[0]
            newdll=dll.split('\\')[-1]
            os.system(f'move {dll} equilipy\\{newdll}')
            newdlls=[newdll]
            self.copy_file(f'equilipy\\{newdll}', os.path.join(build_lib, f'equilipy\\{newdll}'))
            for ext in self.extensions:
                ext.include_package_data=True,
                ext.package_data=dict({'equilipy':newdlls})

        # Continue with the build process
        build_ext.run(self)
        
        # Remove all files
        for f in glob.glob('equilipy/*.mod'):
            os.remove(f)
        for f in glob.glob('equilipy/*.f90'):
            os.remove(f)
        for f in glob.glob('equilipy/*.o'):
            os.remove(f)
        
    def compile_fortran(self, ext):
        # compile
        for file in ext.dirs:
            filename = '/'.join(file.split('/')[1:])
            self.slist=self.slist+f'{filename} '
            
setup(
    name='equilipy',
    version='0.1.1',
    author='Sunyong Kwon',
    author_email='sunyong.kwon@mail.mcgill.ca',
    license='BSD-3',
    packages=find_packages(),
    install_requires=[
        'numpy',
        "glob2",
        "regex",
        "numba",
        "dataclasses",
        "meson",
        "ninja",
        "polars",
        "xlsx2csv",
        "matplotlib",
        "tqdm",
        "mpi4py",
        ],
    ext_modules=[f2py_Extension('equilipy',['equilipy/fsrc'])],  # Include the Fortran extension
    entry_points={
        'console_scripts': [
            'equilipy = equilipy.__main__:main'
        ]
    },
    python_requires='>=3.9',
    cmdclass=dict(build_ext=f2py_Build),  # Use custom build_ext
)
