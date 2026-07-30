#!/usr/bin/env python
import subprocess, sys, os

test_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(test_dir))
c2_target = os.path.join(test_dir, "c2_git")

env = os.environ.copy()

# Ensure c2ImageD11 is available
try:
    import c2ImageD11
    print("c2ImageD11 found in Python path")
except ImportError:
    print("Installing c2ImageD11 from GitHub...")
    subprocess.check_call([
        sys.executable, "-m", "pip", "install",
        "https://github.com/jonwright/c2ImageD11.git",
        "--upgrade",
        "--target={}".format(c2_target),
    ])
    pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = c2_target + (os.pathsep + pythonpath if pythonpath else "")

# Run test suite with f2py backend
print("\n=== Test suite: f2py backend ===")
r1 = subprocess.run(
    [sys.executable, "-m", "pytest", "test/", "-v"],
    cwd=repo_root, env=env,
)

# Run test suite with c2ImageD11 backend
env["IMAGED11_USE_C2"] = "1"
print("\n=== Test suite: c2ImageD11 backend ===")
r2 = subprocess.run(
    [sys.executable, "-m", "pytest", "test/", "-v"],
    cwd=repo_root, env=env,
)

# Summary
print("\n" + "=" * 50)
print("f2py backend:      {} ({})".format(
    'PASS' if r1.returncode == 0 else 'FAIL', r1.returncode))
print("c2ImageD11 backend: {} ({})".format(
    'PASS' if r2.returncode == 0 else 'FAIL', r2.returncode))
sys.exit(0 if r1.returncode == 0 and r2.returncode == 0 else 1)
