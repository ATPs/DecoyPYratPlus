from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
except Exception as exc:
    raise SystemExit('Cython is required to build fast_digest: {}'.format(exc))

extensions = [
    Extension(
        name='fast_digest',
        sources=['decoypyrat/fast_digest.pyx'],
    )
]

setup(
    name='fast_digest',
    ext_modules=cythonize(extensions, compiler_directives={'language_level': '3'}),
)
