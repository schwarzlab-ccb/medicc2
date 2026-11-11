"""ctypes wrapper of the native C++ engine (medicc/native/src).

Each Model owns one native handle with its own caches. Handles are not thread-safe:
every worker thread creates its own Model. Messages (ancestor languages) are immutable
and can be read by other threads while their Model is alive.
"""
import ctypes
import importlib.machinery
import logging
import os
from dataclasses import dataclass

logger = logging.getLogger(__name__)

PARTS_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                          'objects', 'ancestral_parts.txt')
SYMBOLS = set('012345678X')

_library = None


@dataclass(frozen=True)
class ModelConfig:
    """Selects the native model. Returned by ``medicc.io.read_fst()``.

    Allele-specific data (normal CN 1 per allele) and total copy numbers (one column,
    normal CN 2) use the same model; ``total_copy_numbers`` only states the data mode.
    ``no_wgd`` restricts the model to A (losses followed by gains, no WGD).
    The config is picklable, so it can be sent to joblib workers.
    """
    total_copy_numbers: bool = False
    no_wgd: bool = False

    def __post_init__(self):
        for name in ('total_copy_numbers', 'no_wgd'):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")


def library():
    """Load the native engine (built by ``python setup.py build_ext``) once per process."""
    global _library
    if _library is not None:
        return _library
    directory = os.path.dirname(os.path.realpath(__file__))
    paths = [os.path.join(directory, '_engine' + suffix)
             for suffix in importlib.machinery.EXTENSION_SUFFIXES]
    paths = [path for path in paths if os.path.isfile(path)]
    if len(paths) == 0:
        raise RuntimeError("The native MEDICC2 engine is not built. "
                           "Run 'python setup.py build_ext --inplace' or reinstall MEDICC2.")
    lib = ctypes.CDLL(paths[0])
    p, i, s, d = ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_double
    signatures = {
        'am_error': (s, []),
        'nf_load': (p, [s, i]),
        'am_free': (None, [p]),
        'am_string': (p, [s]),
        'am_drop': (None, [p]),
        'am_offset': (d, [p]),
        'am_start': (i, [p]),
        'am_final': (i, [p, i]),
        'am_arc': (i, [p, i, i, ctypes.POINTER(i)]),
        'am_import': (p, [i, i, i, i, ctypes.POINTER(i)]),
        'nf_join': (p, [p, p, p]),
        'nf_resolve': (s, [p, s, p, ctypes.POINTER(d)]),
        'nf_score': (d, [p, s, s, i]),
        'nf_history': (s, [p, s, s]),
        'nf_clear': (None, [p]),
        'nf_finish_row': (None, [p]),
    }
    for name, (result, arguments) in signatures.items():
        function = getattr(lib, name)
        function.restype = result
        function.argtypes = arguments
    _library = lib
    return lib


class Message:
    """Owner of one native ancestor message (a weighted acceptor of CN strings)."""

    def __init__(self, model, handle):
        if not handle:
            raise RuntimeError(model.lib.am_error().decode())
        self.model = model  # keeps the library loaded while the message exists
        self.handle = handle

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
    """One native model handle. Use it as a context manager to release it deterministically.

    Args:
        config (ModelConfig, optional): Model selection. Defaults to the full model.

    Raises:
        RuntimeError: If the engine is not built or the model file cannot be read.
    """

    def __init__(self, config=None):
        self.config = config or ModelConfig()
        if not isinstance(self.config, ModelConfig):
            raise TypeError("config must be a medicc.native.ModelConfig")
        self.handle = None
        self.lib = library()
        self.handle = self.lib.nf_load(PARTS_FILE.encode(), self.config.no_wgd)
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
            raise TypeError("Expected a native message")
        message.check_open()

    def string(self, profile):
        """Return the message that accepts only ``profile`` at cost 0."""
        self.check_open()
        self.check_profile(profile)
        return Message(self, self.lib.am_string(profile.encode()))

    def join(self, left, right):
        """Return the message of all minimum-cost common ancestors of two child messages.

        The message weight of an ancestor is the cost of both subtrees below it.
        """
        self.check_open()
        self.check_message(left)
        self.check_message(right)
        return Message(self, self.lib.nf_join(self.handle, left.handle, right.handle))

    def resolve(self, parent, child):
        """Return one optimal child profile for a fixed parent profile, and its cost.

        Args:
            parent (str): Parent copy-number profile.
            child (Message): Message of the child subtree.

        Returns:
            tuple: (child profile, cost of the branch plus the child subtree)

        Raises:
            RuntimeError: If no child profile is reachable from the parent.
        """
        self.check_open()
        self.check_profile(parent)
        self.check_message(child)
        cost = ctypes.c_double()
        profile = self.lib.nf_resolve(self.handle, parent.encode(), child.handle, ctypes.byref(cost))
        if profile is None:
            raise RuntimeError(self.lib.am_error().decode())
        return profile.decode(), cost.value

    def score(self, x, y, directed=False):
        """Return the native score of two profiles of the same layout.

        Args:
            x (str), y (str): Copy-number profiles, chromosomes separated by X.
            directed (bool): If True, the branch score from parent x to child y;
                otherwise the symmetric distance (x and y descend from a common ancestor).

        Returns:
            float: The score; inf if y cannot be reached.

        Raises:
            ValueError: If the profiles are invalid or their layouts differ.
        """
        self.check_open()
        self.check_profile(x)
        self.check_profile(y)
        if len(x) != len(y) or any((a == 'X') != (b == 'X') for a, b in zip(x, y)):
            raise ValueError("Profiles have different chromosome layouts")
        value = self.lib.nf_score(self.handle, x.encode(), y.encode(), directed)
        if value < 0:
            raise RuntimeError(self.lib.am_error().decode())
        return value

    def history(self, parent, child):
        """Return one minimum-cost ordered event history from ``parent`` to ``child``.

        Returns:
            tuple: (score, stages) where stages is a list of dicts with the keys
            'type' ('loss', 'gain' or 'wgd'), 'before' and 'after' (profiles).

        Raises:
            RuntimeError: If no finite history exists.
        """
        self.check_open()
        self.check_profile(parent)
        self.check_profile(child)
        text = self.lib.nf_history(self.handle, parent.encode(), child.encode())
        if text is None:
            raise RuntimeError(self.lib.am_error().decode())
        lines = text.decode().splitlines()
        return float(lines[0]), [dict(zip(('type', 'before', 'after'), line.split('\t')))
                                 for line in lines[1:]]

    def _finish_row(self):
        """Release row-local FSTs; retain the bounded exact-score cache."""
        self.check_open()
        self.lib.nf_finish_row(self.handle)

    def clear(self):
        """Release the score caches of this handle."""
        self.check_open()
        self.lib.nf_clear(self.handle)

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
