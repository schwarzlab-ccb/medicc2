import os
import pathlib
import subprocess
import time


def test_basic_import():
    "Testing that fstlib can be imported"
    import fstlib


def test_setup():
    "Testing that fstlib can be installed"
    process = subprocess.Popen(["python", "setup.py", "build_ext", "--inplace"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0


def test_setup_deleted_so_cpp_files():
    "Testing that fstlib can be installed after deleting the existing .so and .cpp files"
    delete_process = subprocess.Popen(["rm", "fstlib/cext/*.so", "fstlib/cext/*.cpp"],
                                      stdout=subprocess.PIPE,
                                      cwd=pathlib.Path(__file__).parent.parent.absolute())

    b = pathlib.Path(__file__).parent.parent.absolute()
    a = os.listdir('fstlib/cext')
    c = delete_process.returncode == 0
    
    process = subprocess.Popen(["python", "setup.py", "build_ext", "--inplace"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    c = process.returncode == 0

    assert process.returncode == 0
