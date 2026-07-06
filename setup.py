import glob
import os
from pathlib import Path
import shlex
import subprocess
import sys

from numpy import get_include
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


LAPACK_ARGS_ENV = "EQUILIPY_F2PY_LAPACK_ARGS"
FORTRAN_SOURCE_DIRS = (
    "1_Modules",
    "2_Utils",
    "3_SystemDef",
    "4_GibbsModels",
    "5_Submin",
    "6_ActiveSetLagrangian",
    "7_MultiPhaseMinimizer",
)


def default_lapack_args():
    """Return f2py/Meson link arguments for the platform LAPACK provider."""
    user_args = os.environ.get(LAPACK_ARGS_ENV)
    if user_args:
        return shlex.split(user_args)

    if sys.platform == "darwin":
        veclib = Path(
            "/System/Library/Frameworks/Accelerate.framework"
            "/Versions/Current/Frameworks/vecLib.framework"
        )
        if veclib.exists():
            return [f"-L{veclib}", "-lLAPACK", "-lBLAS"]

    if sys.platform == "win32":
        mingw_prefix = Path(os.environ.get("MINGW_PREFIX", "C:/msys64/mingw64"))
        return [f"-L{mingw_prefix / 'lib'}", "-lopenblas"]

    return ["-llapack", "-lblas"]


def build_environment():
    """Return the environment used by the f2py subprocess."""
    env = os.environ.copy()

    if sys.platform == "win32":
        mingw_prefix = Path(env.get("MINGW_PREFIX", "C:/msys64/mingw64"))
        mingw_bin = str(mingw_prefix / "bin")
        env["PATH"] = f"{mingw_bin}{os.pathsep}{env.get('PATH', '')}"

    return env


class f2py_Extension(Extension):
    """Store Fortran source paths for the custom f2py build."""

    def __init__(self, PackageName, srcdir):
        Extension.__init__(
            self,
            PackageName,
            sources=[],
        )
        sourcedirs = []
        for dirname in FORTRAN_SOURCE_DIRS:
            path = os.path.join(srcdir[0], dirname)
            if not os.path.isdir(path):
                raise FileNotFoundError(f"Missing Fortran source directory: {path}")

            filenames = sorted(glob.glob(f"{path}/*.f*"))
            sourcedirs.extend(filenames)

        self.sourcedirs = [os.path.abspath(sourcedir) for sourcedir in sourcedirs]
        self.dirs = list(sourcedirs)
        self.include_dirs = [get_include()]
        self.packages = find_packages(where=PackageName)


class f2py_Build(build_ext):
    """Build the f2py extension and copy it into the package."""

    slist = ""

    def run(self):
        """Compile Fortran sources and copy the built extension."""
        for ext in self.extensions:
            self.compile_fortran(ext)

        f2py_command = [
            sys.executable,
            "-m",
            "numpy.f2py",
            "-c",
            "-m",
            "equilifort",
            "--backend",
            "meson",
        ]
        f2py_command.extend(self.slist.split())
        f2py_command.extend(default_lapack_args())

        print(f"Running command: {' '.join(f2py_command)}")

        result = subprocess.run(
            f2py_command,
            cwd="equilipy",
            capture_output=True,
            text=True,
            env=build_environment(),
        )

        print("f2py stdout:\n", result.stdout)
        print("f2py stderr:\n", result.stderr)

        if result.returncode != 0:
            raise Exception("f2py compilation failed.")

        lib_name_base = self.get_ext_filename("equilifort").split("/")[-1]

        built_lib_path = None
        for root, _dirs, files in os.walk("equilipy"):
            if lib_name_base in files:
                built_lib_path = os.path.join(root, lib_name_base)
                print(f"Found compiled library at: {built_lib_path}")
                break

        if not built_lib_path:
            raise FileNotFoundError(f"Could not find {lib_name_base} after f2py build.")

        build_lib = os.path.abspath(self.build_lib)
        final_lib_path = os.path.join(build_lib, "equilipy", lib_name_base)

        os.makedirs(os.path.dirname(final_lib_path), exist_ok=True)
        print(f"Copying {built_lib_path} to {final_lib_path}")
        self.copy_file(built_lib_path, final_lib_path)

        self.copy_file(built_lib_path, os.path.join("equilipy", lib_name_base))

        for f in glob.glob("equilipy/*.mod"):
            os.remove(f)

    def compile_fortran(self, ext):
        """Add extension Fortran sources to the f2py command list."""
        for file in ext.dirs:
            rel_path = os.path.relpath(file, "equilipy")
            self.slist += f"{rel_path} "


setup(
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "equilipy.gui": [
            "icons/*.icns",
            "icons/*.ico",
            "icons/*.png",
            "icons/*.svg",
            "fonts/*.otf",
            "fonts/*.ttf",
        ]
    },
    ext_modules=[f2py_Extension("equilipy", ["equilipy/fsrc"])],
    cmdclass=dict(build_ext=f2py_Build),
)
