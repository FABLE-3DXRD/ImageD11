
"""Small compatibility helpers for numpy and its ecosystem (h5py, bslz4_to_sparse)
used across ImageD11."""

import sys
import numpy as np

try:
    import bslz4_to_sparse
except ImportError:
    bslz4_to_sparse = None

# numpy>=2.5 raises a DeprecationWarning on every `a.shape = ...` assignment
# (mutating a possibly-shared array is unsafe, particularly for free-threaded
# builds), regardless of whether the build is actually free-threaded. The
# replacement, reshape(copy=False), only exists from numpy 2.1 -- but that is
# always true here since 2.5 > 2.1.
_SHAPE_SETTER_IS_DEPRECATED = np.lib.NumpyVersion(np.__version__) >= "2.5.0"

if _SHAPE_SETTER_IS_DEPRECATED:
    def reshape_no_copy(a, *shape):
        """Reshape `a` to `shape`, raising ValueError instead of copying."""
        return a.reshape(*shape, copy=False)
else:
    # Plain in-place shape assignment: this is what the codebase always did,
    # raises no warning before numpy 2.5, and needs no minimum numpy version.
    def reshape_no_copy(a, *shape):
        """Reshape `a` to `shape` in place, raising AttributeError instead of copying."""
        a.shape = shape[0] if len(shape) == 1 else shape
        return a


# read_direct_chunk(ds, offset) -> (filter_mask, compressed bytes) of the chunk of
# an h5py dataset starting at offset, e.g. (frame, 0, 0). On python 2, h5py 2.10
# returns the repr of an array.array instead of the bytes. Newer bslz4_to_sparse
# releases read the bytes themselves; without that the repr is parsed back here.
_REPR_HEAD = "array('B', ["


def _h5py_read_direct_chunk(ds, offset):
    return ds.id.read_direct_chunk(offset)


def _py2_read_direct_chunk(ds, offset):
    filter_mask, chunk = ds.id.read_direct_chunk(offset)
    if chunk.startswith(_REPR_HEAD):  # "array('B', [0, 0, 64, ...])"
        chunk = np.fromstring(chunk[len(_REPR_HEAD) : -2], dtype=np.uint8, sep=",").tostring()
    elif chunk == "array('B')":  # empty chunk
        chunk = b""
    return filter_mask, chunk


if bslz4_to_sparse is not None and hasattr(bslz4_to_sparse, "read_direct_chunk"):
    read_direct_chunk = bslz4_to_sparse.read_direct_chunk
elif sys.version_info[0] < 3:
    read_direct_chunk = _py2_read_direct_chunk
else:
    read_direct_chunk = _h5py_read_direct_chunk
