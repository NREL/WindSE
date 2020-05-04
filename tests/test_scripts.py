# import pathlib
# import runpy
# import pytest

# scripts = pathlib.Path(__file__, '..').resolve().glob('*/*.py')

# @pytest.mark.parametrize('script', scripts, ids=lambda script: script.parts[-2]+"/"+script.parts[-1])
# def test_script_execution(script):
#     runpy.run_path(script.as_posix())