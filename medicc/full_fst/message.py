"""Owner of one native DFA pointer; test-only import/export utilities included.

The DFA is allocated by C++ with new. close() calls am_drop exactly once.
Keeping self.model alive also keeps its CDLL loaded while this DFA exists.
"""

import ctypes as C

INF = 100000000


class Message:
    def __init__(self, model, handle):
        self.model = model
        if not handle:
            raise RuntimeError(model.lib.am_error().decode())
        self.handle = handle

    def close(self):
        if self.handle:
            self.model.lib.am_drop(self.handle)
            self.handle = None

    def __del__(self):
        if getattr(self, "handle", None):
            self.close()

    def _open(self):
        if not self.handle:
            raise ValueError("Closed message")

    @property
    def size(self):
        self._open()
        return self.model.lib.am_size(self.handle)

    @property
    def offset(self):
        self._open()
        return self.model.lib.am_offset(self.handle)

    def language(self):
        """Enumerate only for small-case tests; production messages stay compressed."""
        self._open()
        lib = self.model.lib
        result = {}
        stack = [(lib.am_start(self.handle), "", self.offset)]
        while stack:
            s, prefix, cost = stack.pop()
            if s < 0:
                continue
            f = lib.am_final(self.handle, s)
            if f < INF:
                result[prefix] = cost + f
            for c in range(10):
                w = C.c_int()
                t = lib.am_arc(self.handle, s, c, C.byref(w))
                if t >= 0:
                    stack.append(
                        (t, prefix + ("X" if c == 9 else str(c)), cost + w.value)
                    )
        return result

    def to_fst(self, symbols):
        """Export an OpenFST acceptor; put the global offset on final weights."""
        import fstlib

        if not self.handle:
            raise ValueError("Closed message")
        lib = self.model.lib
        f = fstlib.Fst(arc_type="standard")
        f.set_input_symbols(symbols)
        f.set_output_symbols(symbols)
        for _ in range(self.size):
            f.add_state()
        start = lib.am_start(self.handle)
        if start < 0:
            return f
        f.set_start(start)
        for s in range(self.size):
            fin = lib.am_final(self.handle, s)
            if fin < INF:
                f.set_final(s, fstlib.Weight(f.weight_type(), fin + self.offset))
            for c in range(10):
                w = C.c_int()
                t = lib.am_arc(self.handle, s, c, C.byref(w))
                if t >= 0:
                    label = "X" if c == 9 else str(c)
                    f.add_arc(
                        s, (label, label, fstlib.Weight(f.weight_type(), w.value), t)
                    )
        return f

    def choose(self):
        """Copy one minimum-cost profile from the native thread-local output buffer."""
        self._open()
        cost = C.c_int()
        value = self.model.lib.am_choose_profile(self.handle, C.byref(cost))
        if value is None:
            raise RuntimeError(self.model.lib.am_error().decode())
        # c_char_p converts the borrowed C string to Python bytes during the call.
        # decode() gives a Python-owned string; no C pointer escapes this method.
        return value.decode(), cost.value
