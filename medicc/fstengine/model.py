"""ctypes wrapper of the fstengine C++ engine (medicc/fstengine/src).

Each Model owns one engine handle. A handle is not thread-safe: every worker thread
creates its own Model. Close what you open (or use 'with'); __del__ is only a fallback.
"""
import ctypes
import importlib.machinery
import os
from dataclasses import dataclass

# medicc/objects/ancestral_parts.txt, next to this package (written in step 4)
PARTS_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                          'objects', 'ancestral_parts.txt')
SYMBOLS = set('012345678X')

_library = None

@dataclass(frozen=True)
class ModelConfig:
    """Selects the model. Frozen, so it can be hashed and sent to worker processess.

    ''no_wgd'' restricts the model to A (losses followed by gains, no WGD).
    ''total_copy_numbers'' only states the data mode; the model is the same
    """
    total_copy_numbers: bool = False
    no_wgd: bool = False

    def __post_init__(self):
        for name in ('total_copy_numbers', 'no_wgd'):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")

def library():
    """Load _fstengine*.so (next to this file) once per process and declare the C signatures"""
    global _library
    if _library is not None:
        return _library
    directory = os.path.dirname(os.path.realpath(__file__))
    # setup.py names it e.g. _fstengine.cpython-312-x86_64-linux-gnu.so; the hand build _fstengine.so
    paths = [os.path.join(directory, '_fstengine' + suffix)
             for suffix in importlib.machinery.EXTENSION_SUFFIXES]
    paths = [path for path in paths if os.path.isfile(path)]
    if not paths:
        raise RuntimeError("The fstengine is not built. "
                           "Run 'python setup.py build_ext --inplace' or reinstall MEDICC2.")
    lib = ctypes.CDLL(paths[0])
    # Without these declarations ctypes assumes every argument and result is a C int:
    # a 64-bit pointer would be cut to 32 bits, a double read as garbage
    #  short aliases: p = void*, i = int, s = char*, d = double
    p, i, s, d = ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_double

    signatures = { # (return type, [argument type])
        'am_error': (s, []),
        'fstengine_load': (p, [s, i]),
        'am_free': (None, [p]),
        'am_string': (p, [s]),
        'am_drop': (None, [p]),
        'am_offset': (d, [p]),
        'fstengine_join': (p, [p, p, p]),
        'fstengine_resolve': (s, [p, s, p, ctypes.POINTER(d)]),
        'fstengine_score': (d, [p, s, s, i]),
    }
    for name, (result, arguments) in signatures.items():
        function = getattr(lib, name)
        function.restype = result
        function.argtypes = arguments
    _library = lib
    return lib

class Message:
    """Owner of one engine message: the optimal profiles of a tree node, and their cost"""
    """The Python Message doesn't hold the message. It holds a ticket to it. The real Message (fst, offset, length) 
    lives in C++ memory, and Python keeps only its address in self.handle."""

    def __init__(self, model, handle):
        if not handle: # the C function returned nullptr
            raise RuntimeError(model.lib.am_error().decode())
        self.model = model
        self.handle = handle

    @property
    def offset(self):
        self.check_open()
        return self.model.lib.am_offset(self.handle)

    def check_open(self):
        if not self.handle:
            raise ValueError("Closed message")

    def close(self):
        if getattr(self, 'handle', None):
            self.model.lib.am_drop(self.handle)
            self.handle = None

    def __del__(self):
        self.close()


class Model:
    """One engine handle. Use it as a context manager to release it deterministically"""

    def __init__(self, config=None):
        self.config = config or ModelConfig()
        if not isinstance(self.config, ModelConfig):
            raise TypeError("config must be a medicc.fstengine.ModelConfig")
        self.handle = None
        self.lib = library()
        self.handle = self.lib.fstengine_load(PARTS_FILE.encode(), self.config.no_wgd)
        if not self.handle:
            raise RuntimeError(self.lib.am_error().decode())

    def check_open(self):
        if not self.handle:
            raise ValueError("Closed model")

    @staticmethod
    def check_profile(profile):
        if not isinstance(profile, str) or not set(profile) <= SYMBOLS:
            raise ValueError("Expected a string of copy numbers 0..8 and separators X")

    def check_message(self, message):
        if not isinstance(message, Message):
            raise TypeError("Expected an fstengine message")
        message.check_open()

    def string(self, profile):
        """Return the message that accepts only profile at cost 0"""
        self.check_open()
        self.check_profile(profile)
        return Message(self, self.lib.am_string(profile.encode()))

    def join(self, left, right):
        """Return the message of all minimum-cost common ancestors of two child message"""
        self.check_open()
        self.check_message(left)
        self.check_message(right)
        return Message(self, self.lib.fstengine_join(self.handle, left.handle, right.handle))

    def resolve(self, parent, child):
        """Return (one optimal child profile, cost of the branch plus the child subtree).

        Raises RuntimeError if no child profile is reachable from ``parent``.
        """
        self.check_open()
        self.check_profile(parent)
        self.check_message(child)
        cost = ctypes.c_double()
        profile = self.lib.fstengine_resolve(self.handle, parent.encode(), child.handle, ctypes.byref(cost))

        if profile is None:
            raise RuntimeError(self.lib.am_error().decode())
        return profile.decode(), cost.value

    def score(self, x, y, directed=False):
        """T(x, y) if ''directed'' (x is the parent), else K(x,y). inf if impossible"""
        self.check_open()
        self.check_profile(x)
        self.check_profile(y)
        if len(x) != len(y) or any((a == 'X') != (b == 'X') for a, b in zip(x, y)):
            raise ValueError("Profiles have different chromosome layouts.")
        value = self.lib.fstengine_score(self.handle, x.encode(), y.encode(), directed)
        if value < 0:
            raise RuntimeError(self.lib.am_error().decode())
        return value

    def close(self):
        if getattr(self, 'handle', None):
            self.lib.am_free(self.handle)
            self.handle = None

    def __del__(self):
        self.close()

    def __enter__(self):
        self.check_open()
        return self

    def __exit__(self, *exc):
        self.close()









