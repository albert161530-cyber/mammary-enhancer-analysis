from setuptools import Extension, setup
import pybind11

ext_modules = [
    Extension(
        "genoma_kernels",
        ["genoma_kernels.cpp"],
        include_dirs=[pybind11.get_include()],
        extra_compile_args=["/O2"],  # Optimización máxima para MSVC en Windows
    ),
]

setup(
    name="genoma_kernels",
    version="1.0",
    ext_modules=ext_modules,
)