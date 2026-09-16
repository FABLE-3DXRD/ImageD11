
"""Small numpy compatibility helpers used across ImageD11."""

import numpy as np

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
